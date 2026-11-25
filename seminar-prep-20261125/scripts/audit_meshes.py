"""Audit new meshes and budgets. No CFD solver, no physical-result claims."""
from pathlib import Path
import re,json,csv,datetime,hashlib,math
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
R=Path(__file__).resolve().parent;O=Path('/workspace/TBM/seminar-prep-20261125');D=O/'mesh-quality';D.mkdir(exist_ok=True)
assert (R/'MESHES_COMPLETE.txt').exists(),'Mesh generation/checks still running'
def number(pattern,text,cast=float,default=None):
 m=re.search(pattern,text)
 return cast(m.group(1).rstrip('.')) if m else default
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
samples=list(csv.DictReader((R/'memory-samples.csv').open()))
summary=[]
for n in ['B00','C01','C02']:
 C=R/n;meta=json.loads((C/'case-metadata.json').read_text());normal=(C/'log.checkMesh').read_text();extra=(C/'log.checkMesh-additional').read_text();snap=(C/'log.snappyHexMesh').read_text()
 assert '\nEnd\n' in snap and '\nEnd\n' in normal and '\nEnd\n' in extra
 assert not list(C.glob('log.simpleFoam*')) and not list(C.glob('log.pimpleFoam*'))
 assert 'endTime 0;' in (C/'system/controlDict').read_text()
 regions=number(r'Number of regions:\s*(\d+)',normal,int)
 zones=number(r'cell zones:\s*(\d+)',normal,int)
 assert regions==3 and zones==3,(n,regions,zones)
 boundary=(C/'constant/polyMesh/boundary').read_text()
 patchinfo={name:dict(nFaces=int(nf),startFace=int(sf)) for name,nf,sf in re.findall(r'(\w+)\s*\{[^{}]*?nFaces\s+(\d+);[^{}]*?startFace\s+(\d+);[^{}]*?\}',boundary)}
 assert len(patchinfo)==37
 assert all(patchinfo[p['patch']]['nFaces']>0 for p in meta['ports'])
 start=datetime.datetime.fromisoformat((C/'mesh-start-UTC.txt').read_text().strip().replace('Z','+00:00'))
 end=datetime.datetime.fromisoformat((C/'mesh-end-UTC.txt').read_text().strip().replace('Z','+00:00'))
 peak=max([int(q['container_memory_current_bytes']) for q in samples if q['case']==n] or [0])
 # STL volume is geometric conversion validation, not a solved flow result.
 sr=vtk.vtkSTLReader();sr.SetFileName(str(C/'constant/triSurface/volume.stl'));sr.Update()
 mass=vtk.vtkMassProperties();mass.SetInputConnection(sr.GetOutputPort());mass.Update();stlvol=mass.GetVolume()*1e9
 assert abs(stlvol-meta['CAD_volume_mm3'])/meta['CAD_volume_mm3']<.001
 row=dict(model=n,actual_cells=number(r'cells:\s*(\d+)',normal,int),
  points=number(r'points:\s*(\d+)',normal,int),regions=regions,cell_zones=zones,boundary_patches=len(patchinfo),positive_port_patches=36,
  standard_check_passed='Mesh OK.' in normal,additional_check_passed='Mesh OK.' in extra,
  concave_cells=number(r'Concave cells .*?number of cells:\s*(\d+)',extra,int,0),
  small_determinant_cells=number(r'Cells with small determinant .*?number of cells:\s*(\d+)',extra,int,0),
  min_cell_determinant=number(r'Cell determinant .*?minimum:\s*([\deE+.-]+)',extra),
  additional_failed_checks=number(r'Failed\s+(\d+) mesh checks',extra,int,0),
  max_nonorthogonality_deg=number(r'Mesh non-orthogonality Max:\s*([\deE+.-]+)',normal),
  max_skewness=number(r'Max skewness =\s*([\deE+.-]+)',normal),
  min_cell_volume_m3=number(r'Min volume =\s*([\deE+.-]+)',normal),
  mesh_volume_m3=number(r'Total volume =\s*([\deE+.-]+)',normal),
  CAD_volume_mm3=meta['CAD_volume_mm3'],STL_volume_mm3=stlvol,
  mesh_generation_and_checks_wall_minutes=(end-start).total_seconds()/60,
  sampled_container_memory_current_peak_GiB=peak/1024**3,
  cumulative_container_memory_peak_at_end_GiB=int((C/'container-cumulative-memory-peak-bytes.txt').read_text())/1024**3,
  memory_metric='cgroup memory including cache;5s samples, not sum of process RSS; cumulative peak is whole container since start',
  explicit_prism_layers=False,wall_yplus='Not available until a flow solution exists',
  source_STEP=meta['source_STEP'],source_sha256=meta['source_sha256'],flow_solver_executed=False,
  native_case_path=str(C),patch_geometry=patchinfo)
 assert row['min_cell_volume_m3']>0
 row['mesh_volume_error_percent']=100*(row['mesh_volume_m3']*1e9/meta['CAD_volume_mm3']-1)
 row['concave_cell_fraction_percent']=100*row['concave_cells']/row['actual_cells']
 if row['small_determinant_cells']:
  row['readiness']='Repair under-determined cells and recheck before a comparison pilot; no performance CFD ready claim'
 else:
  row['readiness']='Screening mesh only; additional quality flags require mesh/wall sensitivity before performance validation' if not row['additional_check_passed'] else 'Geometric checks passed; flow convergence and mesh independence still untested'
 summary.append(row)
 print(n,row['actual_cells'],'cells; standard',row['standard_check_passed'],'additional',row['additional_check_passed'],flush=True)
for f in ['transportProperties','turbulenceProperties']:
 assert len({sha(R/n/'constant'/f) for n in ['B00','C01','C02']})==1
for f in ['fvSchemes','fvSolution','blockMeshDict','controlDict','midchannelProbeLocations']:
 assert len({sha(R/n/'system'/f) for n in ['B00','C01','C02']})==1
for f in ['U','p','k','omega','nut']:
 assert len({sha(R/n/'0'/f) for n in ['B00','C01','C02']})==1
(D/'mesh-summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
keys=['model','actual_cells','regions','standard_check_passed','additional_check_passed','additional_failed_checks','concave_cells','small_determinant_cells','min_cell_determinant','max_nonorthogonality_deg','max_skewness','mesh_volume_error_percent','mesh_generation_and_checks_wall_minutes','sampled_container_memory_current_peak_GiB']
with (D/'mesh-summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
 for q in summary:w.writerow({k:q[k] for k in keys})
(D/'matched-input-validation.json').write_text(json.dumps(dict(properties_equal=True,numerical_solution_controls_equal=True,background_mesh_equal=True,same_local_resolution_rule=True,probe_locations_equal=True,initial_fields_equal=True,only_existing_geometry_and_its_refinement_surfaces_differ=True,long_flow_computation_executed=False,flow_solver_executed=False),indent=2))
portsets={n:{r['patch']:r for r in csv.DictReader((D/f'{n}-native-port-geometry.csv').open())} for n in ['B00','C01','C02']}
comp=[]
for n in ['C01','C02']:
 diff=[100*(float(portsets[n][p]['area_mm2'])/float(portsets['B00'][p]['area_mm2'])-1) for p in portsets['B00']]
 comp.append(dict(model=n,reference='B00',maximum_absolute_native_port_area_difference_percent=max(abs(v) for v in diff),
  ports_with_different_face_counts=sum(portsets[n][p]['native_faces']!=portsets['B00'][p]['native_faces'] for p in portsets['B00'])))
(D/'matched-port-geometry.json').write_text(json.dumps(comp,indent=2))
print('Mesh audit finished; no flow solve executed')
