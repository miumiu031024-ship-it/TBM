"""Quality-driven internal-face agglomeration, never remove a physical boundary."""
from pathlib import Path
import sys,re,json,shutil,hashlib,argparse
import numpy as np
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/tbm-seminar-prep-20261125')
def read_list(p,dtype):
 b=p.read_bytes();m=re.search(rb'\n\s*(\d+)\s*\(',b);assert m
 return b,np.frombuffer(b,dtype=dtype,offset=m.end(),count=int(m.group(1))),m.end()
parser=argparse.ArgumentParser();parser.add_argument('--tag',choices=['R05'],default='R05');parser.add_argument('models',nargs='*');args=parser.parse_args()
for model in args.models or ['C01']:
 assert model in ['B00','C01','C02']
 old=OLD/model;C=R/f'{model}_{args.tag}';assert not C.exists()
 (C/'constant').mkdir(parents=True);shutil.copytree(old/'constant/polyMesh',C/'constant/polyMesh')
 for f in ['transportProperties','turbulenceProperties']:shutil.copy2(old/'constant'/f,C/'constant'/f)
 shutil.copytree(old/'system',C/'system');shutil.copytree(old/'0',C/'0');(C/'case.foam').touch()
 P=old/'constant/polyMesh';log=(old/'log.checkMesh-additional').read_text()
 m=re.search(r'Cells with small determinant .*?number of cells:\s*(\d+)',log);count=int(m.group(1)) if m else 0
 selected=[];predictions=[]
 if count:
  b=(P/'sets/underdeterminedCells').read_bytes();mm=re.search(rb'\n(\d+)\s*\(([^)]*)\)',b);bad=np.array([int(v) for v in mm.group(2).split()]);assert len(bad)==count
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
  used=set()
  for cid in sorted(map(int,bad)):
   candidates=[];ownfaces=faces_of(cid)
   for fi in ownfaces[ownfaces<len(nei)]:
    other=int(nei[fi] if owner[fi]==cid else owner[fi])
    if other in bad or other in used:continue
    merged=np.setxor1d(ownfaces,faces_of(other));quality=det(merged)
    newcc,pyr=geometry(merged,{cid,other});angle=worst_nonorth(merged,{cid,other},newcc)
    print('candidate',cid,other,'det',quality,'minPyrDistance',min(pyr),'maxNonorth',angle,flush=True)
    if quality>.005 and min(pyr)>1e-12 and angle<85:
     candidates.append((-angle,quality,int(fi),other,float(det(faces_of(other)))))
   assert candidates,cid
   negative_angle,quality,fi,other,neighborq=max(candidates)
   assert quality>.005,(cid,quality)
   selected.append(fi);used.update([cid,other]);predictions.append(dict(original_bad_cell=cid,neighbor_cell=other,
     removed_internal_face=fi,before_bad_determinant=det(ownfaces),before_neighbor_determinant=neighborq,
     predicted_merged_internal_determinant=quality,predicted_max_nonorth=-negative_angle))
  assert len(selected)==count and max(selected)<len(nei)
  text='FoamFile {version 2.0; format ascii; class faceSet; object repairFaces;}\n'+str(len(selected))+'\n(\n'+'\n'.join(map(str,selected))+'\n)\n'
  (C/'constant/polyMesh/sets/repairFaces').write_text(text);(C/'selected-repairFaces.input').write_text(text)
 meta=json.loads((old/'case-metadata.json').read_text());meta.update(revision=args.tag+'_quality_driven_internal_face_cleanup',boundary_merge_cosine=0.999999999999 if args.tag=='R05' else 2,source_case=str(old),
   CAD_geometry_changed=False,flow_solver_executed=False,cleanup_rule='For each checkMesh determinant<.001 cell, merge across one INTERNAL face to a disjoint healthy neighbor with positive face pyramids and minimum predicted maximum nonorthogonality (<85 deg); require determinant>.005. No physical boundary is selected.',
   selected_internal_faces=selected,predictions=predictions)
 (C/'case-metadata.json').write_text(json.dumps(meta,indent=2));print(model,'local cleanup selected',len(selected),'internal faces',predictions,flush=True)
