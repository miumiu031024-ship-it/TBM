"""Reproduce v2412 determinants from native faces, including internal-face mask."""
from pathlib import Path
import re,json,csv,sys
import numpy as np
R=Path(__file__).resolve().parent
OUT=Path('/workspace/TBM/mesh-repair-pilot-20261125');OUT.mkdir(exist_ok=True)
def binary_list(p,dtype,offset=0):
 b=p.read_bytes();m=re.search(rb'\n\s*(\d+)\s*\(',b[offset:]);assert m
 return b,np.frombuffer(b,dtype=dtype,offset=offset+m.end(),count=int(m.group(1))),offset+m.end()
rows=[]
for model in sys.argv[1:] or ['C01','C02']:
 C=R/(model+'_R07');P=C/'constant/polyMesh'
 ss=(P/'sets/repairTargetCells').read_bytes();m=re.search(rb'\n(\d+)\s*\(([^)]*)\)',ss);ids=np.array([int(x) for x in m.group(2).split()]);assert len(ids)==6
 bo,owner,_=binary_list(P/'owner','<i4');bn,nei,_=binary_list(P/'neighbour','<i4')
 bp,flat,_=binary_list(P/'points','<f8');points=flat # read vector-count separately below
 pm=re.search(rb'\n\s*(\d+)\s*\(',bp);points=np.frombuffer(bp,dtype='<f8',offset=pm.end(),count=int(pm.group(1))*3).reshape(-1,3)
 bf,offsets,start=binary_list(P/'faces','<i4');end=start+4*len(offsets)
 cm=re.search(rb'\n\s*(\d+)\s*\(',bf[end:]);assert cm
 conn=np.frombuffer(bf,dtype='<i4',offset=end+cm.end(),count=int(cm.group(1)))
 for cid in ids:
  faces=np.union1d(np.flatnonzero(owner==cid),np.flatnonzero(nei==cid));areas=[];vertices=[]
  for fi in faces:
   idsface=conn[offsets[fi]:offsets[fi+1]];p=points[idsface];vertices.extend(idsface.tolist())
   areas.append(np.cross(p,np.roll(p,-1,axis=0)).sum(axis=0)*.5)
  areas=np.array(areas);mag=np.linalg.norm(areas,axis=1);inside=faces<len(nei);internal=areas[inside]
  avg=np.linalg.norm(internal,axis=1).mean();T=(internal/avg).T@(internal/avg)
  eigen=np.linalg.eigvalsh(T);quality=abs(np.linalg.det(T))/8
  all_tensor=np.einsum('ni,nj,n->ij',areas,areas,1/mag)/mag.sum();snap=np.linalg.det(all_tensor)*27
  pp=points[np.unique(vertices)];centre=pp.mean(0)*1000;s=np.linalg.svd(pp-pp.mean(0),compute_uv=False)
  rows.append(dict(model=model,native_global_cell_id=int(cid),x_mm=float(centre[0]),y_mm=float(centre[1]),z_mm=float(centre[2]),
    faces=len(faces),internal_faces=int(inside.sum()),boundary_faces=int((~inside).sum()),
    checkMesh_internal_face_determinant=float(quality),snappy_all_face_determinant=float(snap),
    internal_tensor_smallest_eigenvalue=float(eigen.min()),internal_tensor_largest_eigenvalue=float(eigen.max()),
    vertex_PCA_min_over_max=float(s.min()/s.max())))
 expected=float(re.search(r'Cell determinant .*?minimum:\s*([\deE+.-]+)',(C/'log.checkMesh-additional').read_text()).group(1))
 got=min(q['checkMesh_internal_face_determinant'] for q in rows if q['model']==model)
 assert abs(got-expected)/expected<1e-4,(model,got,expected)
 print(model,'reproduced minimum determinant',got,'expected',expected,flush=True)
with (OUT/'repaired-target-cell-indicators.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(OUT/'determinant-definition.json').write_text(json.dumps(dict(version='OpenFOAMv2412',
 checkMesh='abs(det(sum outer(Sf/meanInternalArea,Sf/meanInternalArea) over internal/coupled faces))/8',
 snappy='det(sum outer(Sf,Sf)/mag(Sf) over ALL faces /sum mag(Sf))*27',
 source_checkMesh='src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshTools.C:829',
 source_snappy='src/meshTools/primitiveMeshGeometry/primitiveMeshGeometry.C:1127',
 note='Different face selection and normalization; neither is cell volume or a CFD solution residual.'),indent=2))
