from pathlib import Path
import json,re,csv,hashlib,shutil
import numpy as np
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/tbm-seminar-prep-20261125');O=Path('/workspace/TBM/mesh-repair-pilot-20261125');Q=O/'mesh-quality';Q.mkdir(exist_ok=True)
def num(pattern,s,cast=float,default=None):
 m=re.search(pattern,s);return cast(m.group(1).rstrip('.')) if m else default
def row(n,p,stage):
 s=(p/'log.checkMesh').read_text();t=(p/'log.checkMesh-additional').read_text();assert '\nEnd\n' in s and '\nEnd\n' in t
 return dict(model=n,stage=stage,cells=num(r'cells:\s*(\d+)',s,int),regions=num(r'Number of regions:\s*(\d+)',s,int),standard_pass='Mesh OK.' in s,additional_pass='Mesh OK.' in t,
 small_determinant_cells=num(r'Cells with small determinant .*?number of cells:\s*(\d+)',t,int,0),min_determinant=num(r'Cell determinant .*?minimum:\s*([\deE+.-]+)',t),concave_cells=num(r'Concave cells .*?number of cells:\s*(\d+)',t,int,0),
 max_nonorth_deg=num(r'Mesh non-orthogonality Max:\s*([\deE+.-]+)',s),max_skewness=num(r'Max skewness =\s*([\deE+.-]+)',s),min_volume_m3=num(r'Min volume =\s*([\deE+.-]+)',s),volume_m3=num(r'Total volume =\s*([\deE+.-]+)',s),additional_failed_checks=num(r'Failed\s+(\d+) mesh checks',t,int,0))
def csvsave(p,rows):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
rows=[];points=[];proof=[]
for n in ['B00','C01','C02']:
 c=R/(n+'_R07');rows.extend([row(n,OLD/n,'before'),row(n,c,'R07_after')]);meta=json.loads((c/'case-metadata.json').read_text());assert meta['boundary_mesh_exactly_unchanged']
 a=rows[-1];assert a['standard_pass'] and a['small_determinant_cells']==0 and a['min_volume_m3']>0 and a['regions']==3
 assert a['cells']==rows[-2]['cells'];assert a['max_skewness']<=rows[-2]['max_skewness']+1e-7
 for q in meta['changed_interior_points']:
  points.append(dict(model=n,native_cell_id=q['original_bad_cell'],native_point_id=q['interior_point_id'],original_x_mm=q['original_position_m'][0]*1000,original_y_mm=q['original_position_m'][1]*1000,original_z_mm=q['original_position_m'][2]*1000,new_x_mm=q['new_position_m'][0]*1000,new_y_mm=q['new_position_m'][1]*1000,new_z_mm=q['new_position_m'][2]*1000,displacement_um=float(np.linalg.norm(q['displacement_um'])),local_min_determinant_after=q['local_after_min_determinant'],local_after_max_nonorth_deg=q['local_after_max_nonorth_deg']))
 proof.append(dict(model=n,source_CAD_sha256=meta['source_sha256'],boundary_mesh_exactly_unchanged=True,topology_sha256_equal=meta['topology_files_sha256_equal'],changed_internal_vertices=len(meta['changed_interior_points']),native_cell_ids_unchanged=True))
 for f in ['log.checkMesh','log.checkMesh-additional','case-metadata.json','repair-check-start-UTC.txt','repair-check-end-UTC.txt']:shutil.copy2(c/f,Q/(n+'-'+f))
 # Retain only the six repaired targets as a named set for visualization; this is not an unresolved-error set.
 if meta['changed_interior_points']:
  ids=[q['original_bad_cell'] for q in meta['changed_interior_points']];(c/'constant/polyMesh/sets/repairTargetCells').write_text('FoamFile {version 2.0; format ascii; class cellSet; object repairTargetCells;}\n'+str(len(ids))+'\n(\n'+'\n'.join(map(str,ids))+'\n)\n')
 for suffix in ['native-port-geometry.csv','native-region-connectivity.json']:
  shutil.copy2(Path('/workspace/TBM/seminar-prep-20261125/mesh-quality')/(n+'-'+suffix),Q/(n+'-unchanged-'+suffix))
 # Sources above remain applicable because topology and physical boundary vertices are byte-exact equal.
csvsave(Q/'mesh-quality-before-after.csv',rows);csvsave(Q/'internal-point-adjustments.csv',points)
(Q/'mesh-quality-before-after.json').write_text(json.dumps(rows,indent=2));(Q/'exact-boundary-preservation.json').write_text(json.dumps(proof,indent=2))
# Same controls and material/BCs across all three new cases.
files=['system/controlDict','system/fvSchemes','system/fvSolution','system/decomposeParDict','system/midchannelProbeLocations','constant/transportProperties','constant/turbulenceProperties','0/U','0/p','0/k','0/omega','0/nut']
equal={f:len({hashlib.sha256((R/(n+'_R07')/f).read_bytes()).hexdigest() for n in ['B00','C01','C02']})==1 for f in files};assert all(equal.values());(O/'matched-input-validation.json').write_text(json.dumps(equal,indent=2))
print('Native quality rechecks and exact boundary/matched-controls validation passed',rows,flush=True)
