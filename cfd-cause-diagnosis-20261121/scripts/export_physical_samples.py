"""Actual sampled URANS startup animation and point fields, not periodic modes."""
from pathlib import Path
import re,json,os,ast
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
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
R=Path(__file__).resolve().parent;O=R/'diagnosis';rho=746.200955901462
p=R/'analyze_transient.py';tree=ast.parse(p.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='probe_values');exec(compile(ast.Module(body=[node],type_ignores=[]),str(p),'exec'))
all_data={}
for name in ['T30_dt2e5','T40_dt1e5']:
 C=R/name;assert '\nEnd\n' in (C/'log.pimpleFoam').read_text()
 P=C/'postProcessing/interiorProbes/0'
 body=(P/'U').read_text();coords=re.findall(r'^# Probe \d+ \(([^)]+)\)',body,re.M);xyz=np.array([list(map(float,q.split())) for q in coords]);assert xyz.shape==(1359,3)
 pu=probe_values(P/'U',True);pp=probe_values(P/'p');pv=probe_values(P/'vorticity',True);times=sorted(set(pu)&set(pp)&set(pv));assert len(times)==20 and abs(times[-1]-.002)<1e-10
 u=np.stack([pu[t] for t in times]);pressure=np.stack([pp[t] for t in times])*rho;vort=np.stack([pv[t] for t in times]);valid=np.isfinite(u).all(axis=2)&(abs(u)<1e10).all(axis=2)&(abs(vort)<1e10).all(axis=2)
 np.savez_compressed(O/f'{name}-physical-probes.npz',points_m=xyz,time_s=np.array(times),velocity_m_s=u,pressure_Pa=pressure,vorticity_s_inv=vort,valid=valid)
 all_data[name]=(xyz,np.array(times),u,pressure,vort,valid)
 good=valid[-1];pts=vtkPoints();pts.SetData(numpy_to_vtk(xyz[good],deep=True));pd=vtkPolyData();pd.SetPoints(pts)
 for key,val in {'Velocity_m_s':u[-1,good],'Speed_m_s':np.linalg.norm(u[-1,good],axis=1),'Pressure_Pa':pressure[-1,good],'Vorticity_s_inv':vort[-1,good]}.items():
  a=numpy_to_vtk(val,deep=True);a.SetName(key);pd.GetPointData().AddArray(a)
 pd.GetPointData().SetActiveVectors('Velocity_m_s');pd.GetPointData().SetActiveScalars('Speed_m_s')
 g=vtkVertexGlyphFilter();g.SetInputData(pd);g.Update();w=vtkXMLPolyDataWriter();w.SetFileName(str(O/f'{name}-final-probes.vtp'));w.SetInputData(g.GetOutput());w.SetCompressorTypeToZLib();assert w.Write()==1
xyz,times,u,pressure,vort,valid=all_data['T40_dt1e5'];good=valid.all(axis=0)
lon=np.degrees(np.arctan2(xyz[:,1],-xyz[:,0]));lat=np.degrees(np.arcsin(xyz[:,2]/.0464));speed=np.linalg.norm(u,axis=2)
fig,axs=plt.subplots(1,2,figsize=(10,4.6));sc=[]
for ax,vals,cmap,label,lims in zip(axs,[speed[0],pressure[0]/1000],['turbo','viridis'],['Speed [m/s]','Gauge static pressure [kPa]'],[(0,1.5),(float(pressure[:,good].min()/1000),float(pressure[:,good].max()/1000))]):
 a=ax.scatter(lon[good],lat[good],c=vals[good],s=10,cmap=cmap,vmin=lims[0],vmax=lims[1],linewidth=0);sc.append(a);fig.colorbar(a,ax=ax,shrink=.7,label=label);ax.set(xlabel='Longitude [deg]',ylabel='Elevation [deg]',xlim=(-65,65),ylim=(0,65))
title=fig.suptitle('')
fig.subplots_adjust(top=.78,wspace=.30)
def frame(i):
 sc[0].set_array(speed[i,good]);sc[1].set_array(pressure[i,good]/1000)
 title.set_text(f'Actual URANS startup probes, physical t={times[i]:.5f} s\nLeft region r=46.4 mm, dt=1e-5 s; NOT a periodic cycle or converged mean')
 return [*sc,title]
anim=FuncAnimation(fig,frame,frames=len(times),interval=150,blit=False);anim.save(O/'physical-startup-probes.gif',writer=PillowWriter(fps=6),dpi=110);frame(len(times)-1);fig.savefig(O/'physical-final-probes.png',dpi=170);plt.close(fig)
print('Exported actual 20-frame, 1359-point histories, GIF and VTP files')
