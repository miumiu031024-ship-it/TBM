from pathlib import Path
import sys,shutil,json,re,hashlib,numpy as np
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/tbm-seminar-prep-20261125')
for n in sys.argv[1:]:
 src=OLD/n;c=R/(n+'_R07');assert not c.exists();(c/'constant').mkdir(parents=True)
 shutil.copytree(src/'constant/polyMesh',c/'constant/polyMesh')
 for f in ['transportProperties','turbulenceProperties']:shutil.copy2(src/'constant'/f,c/'constant'/f)
 shutil.copytree(src/'system',c/'system');shutil.copytree(src/'0',c/'0');(c/'case.foam').touch()
 adjustments=[] if n=='B00' else json.loads((R/(n+'-R07-point-adjustments.json')).read_text())
 p=c/'constant/polyMesh/points';b=p.read_bytes();m=re.search(rb'\n\s*(\d+)\s*\(',b);a=np.frombuffer(b,dtype='<f8',offset=m.end(),count=3*int(m.group(1))).reshape(-1,3).copy();before=a.copy()
 for q in adjustments:
  i=q['interior_point_id'];assert np.array_equal(a[i],q['original_position_m']);a[i]=q['new_position_m']
 changed=np.flatnonzero(np.any(a!=before,axis=1));assert len(changed)==len(adjustments)
 if adjustments:p.write_bytes(b[:m.end()]+a.tobytes()+b[m.end()+a.nbytes:])
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 equal={f:sha(src/'constant/polyMesh'/f)==sha(c/'constant/polyMesh'/f) for f in ['faces','owner','neighbour','boundary','cellZones','faceZones','pointZones']};assert all(equal.values())
 # A changed vertex is prohibited from belonging to ANY physical boundary face.
 bf=(c/'constant/polyMesh/faces').read_bytes();fm=re.search(rb'\n\s*(\d+)\s*\(',bf);offs=np.frombuffer(bf,dtype='<i4',offset=fm.end(),count=int(fm.group(1)));end=fm.end()+offs.nbytes;cm=re.search(rb'\n\s*(\d+)\s*\(',bf[end:]);conn=np.frombuffer(bf,dtype='<i4',offset=end+cm.end(),count=int(cm.group(1)))
 nb=(c/'constant/polyMesh/neighbour').read_bytes();nm=re.search(rb'\n\s*(\d+)\s*\(',nb);ni=int(nm.group(1));assert not np.isin(changed,conn[offs[ni]:]).any()
 meta=json.loads((src/'case-metadata.json').read_text());meta.update(revision='R07_internal_vertex_quality_optimization',source_case=str(src),CAD_geometry_changed=False,boundary_mesh_exactly_unchanged=True,topology_files_sha256_equal=equal,changed_interior_points=adjustments,flow_solver_executed=False,repair_rule='Only one strictly internal vertex per determinant<.001 cell. Minimize displacement, local determinants>=.0011, positive pyramids, affected internal nonorthogonality<=65deg. No boundary vertices moved.')
 (c/'case-metadata.json').write_text(json.dumps(meta,indent=2));print(n,'R07 copied; moved',len(changed),'interior points; exact physical boundary preservation verified',flush=True)
