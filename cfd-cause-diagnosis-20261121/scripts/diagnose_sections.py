"""Actual baseline field: pressure/velocity at outlet pipe sections, no flow claim from colour."""
from pathlib import Path
import json,csv,os,sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkCommonDataModel import vtkPlane
from vtkmodules.vtkFiltersCore import vtkCutter,vtkTriangleFilter
from vtkmodules.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;O=R/'diagnosis';O.mkdir(exist_ok=True)
M=json.loads((R/'SeedSide/case-metadata.json').read_text());rho=M['rho_kg_m3']
name=sys.argv[1] if len(sys.argv)>1 else 'SeedSide'
series=json.loads((R/name/'VTK'/f'{name}.vtm.series').read_text())['files']
s=max(series,key=lambda s:s['time']); snapshot=s['time']; d=R/name/'VTK'/Path(s['name']).stem/'internal.vtu'
r=vtkXMLUnstructuredGridReader();r.SetFileName(str(d));r.Update();data=r.GetOutput();assert data.GetNumberOfCells()==1264925
ports=[p for p in M['ports'] if p['region']=='Water_Xminus' and p['side']=='OUT']
rows=[]
for zmm in [-9,-5,-1]:
 plane=vtkPlane();plane.SetOrigin(0,0,zmm*.001);plane.SetNormal(0,0,1)
 cut=vtkCutter();cut.SetCutFunction(plane);cut.SetInputData(data)
 tri=vtkTriangleFilter();tri.SetInputConnection(cut.GetOutputPort());tri.Update();d=tri.GetOutput()
 t=vtk_to_numpy(d.GetPolys().GetConnectivityArray()).reshape(-1,3);v=vtk_to_numpy(d.GetPoints().GetData())[t]
 area=np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),axis=1)/2
 xyz=v.mean(axis=1);u=vtk_to_numpy(d.GetCellData().GetArray('U'));p=vtk_to_numpy(d.GetCellData().GetArray('p'))*rho
 for port in ports:
  mask=np.linalg.norm(xyz[:,:2]-np.array([port['x_mm'],port['y_mm']])*.001,axis=1)<.00155
  a=area[mask];uv=u[mask];pa=p[mask];un=-uv[:,2];q=a*un*rho
  assert len(a)>100 and abs(a.sum()/7.068583470577e-6-1)<.02,(port,zmm,a.sum())
  totalp=pa+.5*rho*(uv**2).sum(axis=1)
  rows.append(dict(outlet=port['patch'],z_mm=zmm,section_area_mm2=float(a.sum()*1e6),area_mean_static_Pa=float(np.average(pa,weights=a)),outward_mass_kg_s=float(q.sum()),reverse_mass_kg_s=float(-q[q<0].sum()),area_mean_axial_m_s=float(np.average(un,weights=a)),transverse_velocity_RMS_m_s=float(np.sqrt(np.average((uv[:,:2]**2).sum(axis=1),weights=a))),reverse_area_fraction=float(a[un<0].sum()/a.sum()),flux_weighted_total_Pa=float(np.dot(totalp,q)/q.sum())))
with (O/(name+'-pipe-sections.csv')).open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fig,axs=plt.subplots(1,3,figsize=(12,3.8))
for port in ports:
 q=[v for v in rows if v['outlet']==port['patch']]
 for ax,col in zip(axs,['area_mean_static_Pa','outward_mass_kg_s','flux_weighted_total_Pa']):
  ax.plot([t['z_mm'] for t in q],[t[col]*(1000 if 'kg_s' in col else .001) for t in q],'o-',label=port['patch'])
axs[0].set_ylabel('Area mean static pressure [kPa]');axs[1].set_ylabel('Section mass flux [g/s]');axs[2].set_ylabel('Flux-weighted total pressure [kPa]')
for ax in axs:ax.set_xlabel('Z in outlet pipe [mm]');ax.grid(alpha=.2)
axs[-1].legend(fontsize=7)
fig.suptitle(f'{name}, snapshot {snapshot:g}: actual outlet-pipe sections\nPilot snapshot; sectional U-derived flux is not native face flux phi',fontsize=11)
fig.tight_layout();fig.savefig(O/(name+'-pipe-pressure-flow.png'),dpi=180);plt.close(fig)
print(json.dumps(rows,indent=2))
