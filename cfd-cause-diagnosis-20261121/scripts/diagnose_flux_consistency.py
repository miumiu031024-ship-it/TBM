from pathlib import Path
import ast,re,json,csv
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.vtkFiltersCore import vtkTriangleFilter
from vtkmodules.util.numpy_support import vtk_to_numpy,numpy_to_vtk
R=Path(__file__).resolve().parent;rho=746.200955901462
# Reuse the previously verified binary OpenFOAM patch reader without importing
# the analysis program's top-level execution.
p=Path('/workspace/scratch/tbm-A01-cfd-20261120/analyze.py');tree=ast.parse(p.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='scalar_patch');exec(compile(ast.Module(body=[node],type_ignores=[]),str(p),'exec'))
rows=[]
for name in ['SeedSide',*json.loads((R/'trials.json').read_text()),'T30_dt2e5','T40_dt1e5']:
 C=R/name;series=C/'VTK'/f'{name}.vtm.series'
 if not series.exists():continue
 entries=json.loads(series.read_text())['files'];entry=max(entries,key=lambda x:x['time']);time=entry['time'];folder=C/'VTK'/Path(entry['name']).stem
 if not (C/f'{time:g}'/'phi').exists():continue
 for i in range(1,7):
  patch=f'outlet{i:02d}';r=vtkXMLPolyDataReader();r.SetFileName(str(folder/'boundary'/f'{patch}.vtp'));r.Update();d=r.GetOutput();nf=d.GetNumberOfCells()
  ids=numpy_to_vtk(np.arange(nf),deep=True);ids.SetName('nativeFaceId');d.GetCellData().AddArray(ids)
  tf=vtkTriangleFilter();tf.SetInputData(d);tf.Update();t=tf.GetOutput();xyz=vtk_to_numpy(t.GetPoints().GetData());ff=vtk_to_numpy(t.GetPolys().GetConnectivityArray()).reshape(-1,3)
  av=np.cross(xyz[ff[:,1]]-xyz[ff[:,0]],xyz[ff[:,2]]-xyz[ff[:,0]])*.5; assert (av[:,2]<0).all()
  fi=vtk_to_numpy(t.GetCellData().GetArray('nativeFaceId')).astype(int);A=np.column_stack([np.bincount(fi,weights=av[:,k],minlength=nf) for k in range(3)])
  u=vtk_to_numpy(d.GetCellData().GetArray('U'));qu=rho*np.sum(u*A,axis=1);qp=scalar_patch(C/f'{time:g}'/'phi',patch,nf)*rho
  pressure=vtk_to_numpy(d.GetCellData().GetArray('p'))
  assert np.max(abs(pressure))<1e-12,(name,patch,'Outlet gauge pressure is not zero')
  rows.append(dict(case=name,snapshot=time,outlet=patch,native_phi_kg_s=float(qp.sum()),from_boundary_U_kg_s=float(qu.sum()),difference_kg_s=float(qp.sum()-qu.sum()),native_face_reverse_kg_s=float(-qp[qp<0].sum()),boundary_U_reverse_kg_s=float(-qu[qu<0].sum())))
with (R/'diagnosis/exit-flux-consistency.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
