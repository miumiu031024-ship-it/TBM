from pathlib import Path
import numpy as np,re,json,itertools,sys
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/tbm-seminar-prep-20261125')
def read_list(p,dtype):
 b=p.read_bytes();m=re.search(rb'\n\s*(\d+)\s*\(',b)
 return b,np.frombuffer(b,dtype=dtype,offset=m.end(),count=int(m.group(1))),m.end()
model=sys.argv[1];P=OLD/model/'constant/polyMesh'
b=(P/'sets/underdeterminedCells').read_bytes();mm=re.search(rb'\n(\d+)\s*\(([^)]*)\)',b);bad=np.array([int(v) for v in mm.group(2).split()]);assert len(bad)==6
bo,owner,_=read_list(P/'owner','<i4');bn,nei,_=read_list(P/'neighbour','<i4')
bp,_,start=read_list(P/'points','<f8');pm=re.search(rb'\n\s*(\d+)\s*\(',bp);points=np.frombuffer(bp,dtype='<f8',offset=pm.end(),count=3*int(pm.group(1))).reshape(-1,3).copy()
bf,offsets,start=read_list(P/'faces','<i4');end=start+4*len(offsets);fm=re.search(rb'\n\s*(\d+)\s*\(',bf[end:]);assert fm
conn=np.frombuffer(bf,dtype='<i4',offset=end+fm.end(),count=int(fm.group(1)))
facecache={};areacache={};centercache={};cellcache={}
def faces_of(cid):
 if cid not in facecache:facecache[cid]=np.union1d(np.flatnonzero(owner==cid),np.flatnonzero(nei==cid))
 return facecache[cid]
def area(fi):
 fi=int(fi)
 if fi not in areacache:
  p=points[conn[offsets[fi]:offsets[fi+1]]];p=p-p.mean(0);areacache[fi]=np.cross(p,np.roll(p,-1,axis=0)).sum(0)*.5
 return areacache[fi]
def fc(fi):
 fi=int(fi)
 if fi not in centercache:
  p=points[conn[offsets[fi]:offsets[fi+1]]];mid=p.mean(0);w=np.linalg.norm(np.cross(p-mid,np.roll(p,-1,axis=0)-mid),axis=1)
  centercache[fi]=np.sum(w[:,None]*(p+np.roll(p,-1,axis=0)+mid)/3,axis=0)/sum(w)
 return centercache[fi]
def geometry(faces,pair):
 centers=np.array([fc(f) for f in faces]);A=np.array([area(f)*(1 if int(owner[f]) in pair else -1) for f in faces]);est=centers.mean(0)
 v=np.sum(A*(centers-est),axis=1);cc=np.sum(v[:,None]*(.75*centers+.25*est),axis=0)/sum(v)
 pyr=np.sum(A*(centers-cc),axis=1)/np.linalg.norm(A,axis=1)
 return cc,pyr
def cc_of(cid):
 if cid not in cellcache:cellcache[cid]=geometry(faces_of(cid),{cid})[0]
 return cellcache[cid]
def worst_nonorth(faces,pair,cc):
 values=[]
 for f in faces:
  if f>=len(nei):continue
  other=int(nei[f] if int(owner[f]) in pair else owner[f]);d=cc_of(other)-cc;a=area(f)*(1 if int(owner[f]) in pair else -1)
  values.append(float(np.degrees(np.arccos(np.clip(np.dot(d,a)/np.linalg.norm(d)/np.linalg.norm(a),-1,1)))))
 return max(values)
def det(faces):
 A=np.array([area(f) for f in faces if f<len(nei)]);avg=np.linalg.norm(A,axis=1).mean();T=(A/avg).T@(A/avg)
 return abs(float(np.linalg.det(T)))/8


from scipy.optimize import minimize
boundaryvertices=set(map(int,np.unique(conn[offsets[len(nei)]:])))
records=[]
for cid in sorted(map(int,bad)):
 fs=faces_of(cid);vs=np.unique(np.concatenate([conn[offsets[f]:offsets[f+1]] for f in fs]));interior=[int(v) for v in vs if int(v) not in boundaryvertices];assert len(interior)==1
 v=interior[0];original=points[v].copy()
 positions=np.flatnonzero(conn==v);incidentfaces=np.unique(np.searchsorted(offsets,positions,side='right')-1);assert max(incidentfaces)<len(nei)
 cells=sorted(set(map(int,owner[incidentfaces]))|set(map(int,nei[incidentfaces])))
 allfaces=np.unique(np.concatenate([faces_of(c) for c in cells]));intfaces=allfaces[allfaces<len(nei)]
 # Warm immutable neighboring geometry cache before coordinate changes.
 neighbors=sorted(set(map(int,owner[intfaces]))|set(map(int,nei[intfaces])))
 for c in neighbors:cc_of(c)
 scale=1e-6
 def assess(x):
  points[v]=original+np.array(x)*scale
  for f in incidentfaces:areacache.pop(int(f),None);centercache.pop(int(f),None)
  for c in cells:cellcache.pop(c,None)
  qs=np.array([det(faces_of(c)) for c in cells]);pyrs=[]
  for c in cells:
   cc,pyr=geometry(faces_of(c),{c});cellcache[c]=cc;pyrs.extend(pyr)
  cosines=[]
  for f in intfaces:
   d=cc_of(int(nei[f]))-cc_of(int(owner[f]));a=area(f)
   cosines.append(np.dot(d,a)/np.linalg.norm(d)/np.linalg.norm(a))
  return qs,np.array(pyrs),np.array(cosines)
 base=assess([0,0,0]);print('BEFORE',cid,v,'q',min(base[0]),'pyr',min(base[1]),'angle',np.degrees(np.arccos(min(base[2]))),flush=True)
 def constraints(x):
  q,p,c=assess(x)
  return np.r_[q/.001-1.1,p/scale-.01,(c-np.cos(np.radians(65)))*10]
 opt=minimize(lambda x:float(np.dot(x,x)),np.zeros(3),method='SLSQP',bounds=[(-25,25)]*3,constraints=[dict(type='ineq',fun=constraints)],options=dict(maxiter=180,ftol=1e-10,eps=1e-4))
 q,p,c=assess(opt.x);valid=min(constraints(opt.x))> -1e-5
 print('AFTER',cid,'success',opt.success,'valid',valid,'delta_um',opt.x,'q',min(q),'pyr',min(p),'angle',np.degrees(np.arccos(min(c))),'message',opt.message,flush=True)
 assert valid,(cid,opt)
 records.append(dict(original_bad_cell=cid,interior_point_id=v,original_position_m=original.tolist(),new_position_m=points[v].tolist(),displacement_um=opt.x.tolist(),affected_cells=cells,local_before_min_determinant=float(min(base[0])),local_after_min_determinant=float(min(q)),local_after_max_nonorth_deg=float(np.degrees(np.arccos(min(c)))),optimizer_success=bool(opt.success)))
(R/(model+'-R07-point-adjustments.json')).write_text(json.dumps(records,indent=2))
print('SUCCESS',len(records),'interior vertices; boundary vertices unchanged',flush=True)
