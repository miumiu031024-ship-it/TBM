"""Visualize actual steady-iteration fields, never as physical-time evolution."""
from pathlib import Path
import json,re,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-ramp-plots')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-ramp-cache')
Path(os.environ['XDG_CACHE_HOME']).mkdir(parents=True,exist_ok=True)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkFiltersGeneral import vtkVertexGlyphFilter
from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
from vtkmodules.util.numpy_support import numpy_to_vtk
R=Path(__file__).resolve().parent;O=R/'results';rho=746.200955901462
xyz=None;records=[];history=[]
for name in ['F01_upwind_long','F02_second_order_restore']:
 values={field:{} for field in ['U','p']}
 for field in values:
  for p in sorted((R/name/'postProcessing/interiorProbes').glob('*/'+field),key=lambda p:float(p.parent.name)):
   body=p.read_text()
   if xyz is None:xyz=np.array([list(map(float,m.split())) for m in re.findall(r'^# Probe \d+ \(([^)]+)\)',body,re.M)]);assert xyz.shape==(1359,3)
   for line in body.splitlines():
    if not line.strip() or line.startswith('#'):continue
    a=np.fromstring(line.replace('(',' ').replace(')',' '),sep=' ')
    count=1+1359*(3 if field=='U' else 1)
    if len(a)==count:values[field][a[0]]=a[1:].reshape(-1,3) if field=='U' else a[1:]
 for t in sorted(set(values['U'])&set(values['p'])):
  u=values['U'][t];p=values['p'][t]*rho;assert np.isfinite(u).all() and (abs(u)<1e10).all()
  records.append((name,t,u,p));history.append(dict(case=name,iteration=t))
 U=np.stack([q[2] for q in records if q[0]==name]);P=np.stack([q[3] for q in records if q[0]==name]);T=np.array([q[1] for q in records if q[0]==name]);np.savez_compressed(O/f'{name}-monitor-fields.npz',points_m=xyz,iterations=T,velocity_m_s=U,pressure_Pa=P,physical_time=False)
lon=np.degrees(np.arctan2(xyz[:,1],-xyz[:,0]));lat=np.degrees(np.arcsin(xyz[:,2]/.0464));speed=np.stack([np.linalg.norm(q[2],axis=1) for q in records]);pressure=np.stack([q[3] for q in records])/1000
fig,axs=plt.subplots(1,2,figsize=(10,4.6));dots=[]
for ax,vals,cmap,clim,label in zip(axs,[speed[0],pressure[0]],['turbo','viridis'],[(0,1.5),(float(pressure.min()),float(pressure.max()))],['Speed [m/s]','Gauge static pressure [kPa]']):
 a=ax.scatter(lon,lat,c=vals,s=10,cmap=cmap,vmin=clim[0],vmax=clim[1],linewidth=0);dots.append(a);fig.colorbar(a,ax=ax,shrink=.7,label=label);ax.set(xlabel='Longitude [deg]',ylabel='Elevation [deg]',xlim=(-65,65),ylim=(0,65))
title=fig.suptitle('');fig.subplots_adjust(top=.78,wspace=.3)
def frame(i):
 name,t,u,p=records[i];dots[0].set_array(speed[i]);dots[1].set_array(pressure[i]);title.set_text(f'Actual steady-iteration probes: {name}, iteration {t:g}\nNumerical iteration sequence only; NOT physical time or a vortex period');return [*dots,title]
an=FuncAnimation(fig,frame,frames=len(records),interval=120,blit=False);an.save(O/'numerical-iteration-history.gif',writer=PillowWriter(fps=8),dpi=105);plt.close(fig)
(O/'monitor-frame-index.json').write_text(json.dumps(dict(points=1359,frames=history,physical_time=False),indent=2))
for name in ['F01_upwind_long','F02_second_order_restore']:
 d=np.load(O/f'{name}-r46.4.npz');v=d['valid'];pts=vtkPoints();pts.SetData(numpy_to_vtk(d['points_m'][v],deep=True));pd=vtkPolyData();pd.SetPoints(pts)
 for key,val in {'Velocity_m_s':d['velocity_m_s'][v],'Speed_m_s':np.linalg.norm(d['velocity_m_s'][v],axis=1),'Pressure_Pa':d['pressure_Pa'][v]}.items():
  a=numpy_to_vtk(val,deep=True);a.SetName(key);pd.GetPointData().AddArray(a)
 pd.GetPointData().SetActiveVectors('Velocity_m_s');pd.GetPointData().SetActiveScalars('Speed_m_s');g=vtkVertexGlyphFilter();g.SetInputData(pd);g.Update();w=vtkXMLPolyDataWriter();w.SetFileName(str(O/f'{name}-final-interior-points.vtp'));w.SetInputData(g.GetOutput());w.SetCompressorTypeToZLib();assert w.Write()==1
print('Exported actual numerical-iteration history and final sampled fields')
