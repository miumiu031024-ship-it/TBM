from pathlib import Path
import numpy as np,re,json,itertools,sys
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/tbm-seminar-prep-20261125')
def read_list(p,dtype):
 b=p.read_bytes();m=re.search(rb'\n\s*(\d+)\s*\(',b)
 return b,np.frombuffer(b,dtype=dtype,offset=m.end(),count=int(m.group(1))),m.end()
model=sys.argv[1];P=OLD/model/'constant/polyMesh'
b=(P/'sets/underdeterminedCells').read_bytes();mm=re.search(rb'\n(\d+)\s*\(([^)]*)\)',b);bad=np.array([int(v) for v in mm.group(2).split()]);assert len(bad)==6
bo,owner,_=read_list(P/'owner','<i4');bn,nei,_=read_list(P/'neighbour','<i4')
bp,_,start=read_list(P/'points','<f8');pm=re.search(rb'\n\s*(\d+)\s*\(',bp);points=np.frombuffer(bp,dtype='<f8',offset=pm.end(),count=3*int(pm.group(1))).reshape(-1,3)
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

used=set();clusters=[]
for cid in sorted(map(int,bad)):
 beam=[frozenset([cid])];visited=set(beam);candidates=[]
 for size in range(2,9):
  scored=[]
  for orig in beam:
   of=np.array([],dtype=int)
   for i in orig:of=np.setxor1d(of,faces_of(i))
   adjacent={int(nei[f] if int(owner[f]) in orig else owner[f]) for f in of if f<len(nei)}
   for other in adjacent:
    if other in used or other in bad:continue
    pair=orig|{other}
    if pair in visited:continue
    visited.add(pair);fs=np.setxor1d(of,faces_of(other));q=det(fs);newcc,pyr=geometry(fs,pair);angle=worst_nonorth(fs,pair,newcc)
    score=(min(pyr)*1e6- max(0,angle-65),q)
    scored.append((score,pair))
    if q>.005 and min(pyr)>1e-12 and angle<68:
     candidates.append((-angle,q,sorted(pair),newcc.tolist()))
  print('search',cid,'size',size,'tested',len(scored),'valid',len(candidates),'bestscore',max([x[0] for x in scored],default=None),flush=True)
  if candidates:break
  beam=[x[1] for x in sorted(scored,key=lambda x:x[0],reverse=True)[:6]]
 if not candidates:raise RuntimeError('No safe cluster '+str(cid))
 best=max(candidates);pair=set(best[2]);used.update(pair)
 removed=[int(f) for i in pair for f in faces_of(i) if f<len(nei) and int(owner[f]) in pair and int(nei[f]) in pair]
 clusters.append(dict(original_bad_cell=cid,cells=sorted(pair),removed_faces=sorted(set(removed)),predicted_max_nonorth=-best[0],predicted_merged_determinant=best[1],new_center_m=best[3]))
(R/(model+'-R06-clusters.json')).write_text(json.dumps(clusters,indent=2));print('SELECTED',clusters,flush=True)
