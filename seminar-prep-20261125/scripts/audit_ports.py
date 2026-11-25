"""Pure native boundary geometry audit. Field arrays disabled; no solver."""
from pathlib import Path
import sys,json,csv,math
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
R=Path(__file__).resolve().parent;D=Path('/workspace/TBM/seminar-prep-20261125/mesh-quality');D.mkdir(exist_ok=True)
def visit(data):
 if isinstance(data,vtk.vtkMultiBlockDataSet):
  for i in range(data.GetNumberOfBlocks()):
   b=data.GetBlock(i)
   if b is None:continue
   name=data.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
   if isinstance(b,vtk.vtkPolyData):yield name,b
   else:yield from visit(b)
for n in sys.argv[1:] or ['B00','C01','C02']:
 C=R/n;assert (C/'mesh-end-UTC.txt').exists()
 meta=json.loads((C/'case-metadata.json').read_text());ports=meta['ports']
 reader=vtk.vtkOpenFOAMReader();reader.SetFileName(str(C/'case.foam'));reader.SetCreateCellToPoint(False);reader.SetReadZones(False);reader.SetSkipZeroTime(False);reader.UpdateInformation()
 for i in range(reader.GetNumberOfCellArrays()):reader.SetCellArrayStatus(reader.GetCellArrayName(i),0)
 for i in range(reader.GetNumberOfPatchArrays()):
  name=reader.GetPatchArrayName(i);reader.SetPatchArrayStatus(name,int(name.startswith('patch/inlet') or name.startswith('patch/outlet')))
 reader.Update();patches=dict(visit(reader.GetOutput()));rows=[]
 for p in ports:
  patch=p['patch'];data=patches[patch];assert data.GetNumberOfCells()>0
  assert data.GetCellData().GetNumberOfArrays()==0,'Unexpected field arrays; this audit should be geometry only'
  tri=vtk.vtkTriangleFilter();tri.SetInputData(data);tri.Update();td=tri.GetOutput()
  points=vtk_to_numpy(td.GetPoints().GetData()).astype(float)
  ids=vtk_to_numpy(td.GetPolys().GetConnectivityArray()).reshape(-1,3);vv=points[ids]
  areas=np.linalg.norm(np.cross(vv[:,1]-vv[:,0],vv[:,2]-vv[:,0]),axis=1)/2
  area=float(areas.sum());center=np.average(vv.mean(axis=1),axis=0,weights=areas)
  expected=math.pi*.0015**2;error=100*(area/expected-1)
  centreerr=np.linalg.norm(center-np.array([p[k] for k in ['x_mm','y_mm','z_mm']])*.001)*1000
  row=dict(model=n,patch=patch,CAD_port=p['port'],region=p['region'],native_faces=data.GetNumberOfCells(),area_mm2=area*1e6,CAD_area_mm2=expected*1e6,area_error_percent=error,centroid_error_mm=float(centreerr))
  rows.append(row)
  if patch=='outlet06':
   pts=vtk_to_numpy(data.GetPoints().GetData()).astype(float)
   polygons=[]
   for j in range(data.GetNumberOfCells()):
    cell=data.GetCell(j);polygons.append([[float(pts[cell.GetPointId(k),i]*1000-p[['x_mm','y_mm'][i]]) for i in [0,1]] for k in range(cell.GetNumberOfPoints())])
   (D/f'{n}-OUT06-native-mesh.json').write_text(json.dumps(dict(model=n,native_faces=data.GetNumberOfCells(),polygons_local_mm=polygons,not_CFD=True),indent=2))
 with (D/f'{n}-native-port-geometry.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 assert len(rows)==36
 checks=dict(model=n,ports=36,all_areas_positive=all(q['area_mm2']>0 for q in rows),
  maximum_absolute_area_error_percent=max(abs(q['area_error_percent']) for q in rows),
  maximum_centroid_error_mm=max(q['centroid_error_mm'] for q in rows),
  selected_field_arrays=0,flow_solver_executed=False)
 (D/f'{n}-port-geometry-summary.json').write_text(json.dumps(checks,indent=2))
 print(checks,flush=True)
