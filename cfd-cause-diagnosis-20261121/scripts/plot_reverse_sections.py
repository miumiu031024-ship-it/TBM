"""Actual signed velocity in internal outlet sections, not exit boundary backflow."""
from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkCommonDataModel import vtkPlane
from vtkmodules.vtkFiltersCore import vtkCutter,vtkTriangleFilter
from vtkmodules.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
R=Path(__file__).resolve().parent;M=json.loads((R/'SeedSide/case-metadata.json').read_text())
r=vtkXMLUnstructuredGridReader();r.SetFileName(str(R/'SeedSide/VTK/SeedSide_0/internal.vtu'));r.Update();data=r.GetOutput()
plane=vtkPlane();plane.SetOrigin(0,0,-.001);plane.SetNormal(0,0,1)
cut=vtkCutter();cut.SetCutFunction(plane);cut.SetInputData(data)
tri=vtkTriangleFilter();tri.SetInputConnection(cut.GetOutputPort());tri.Update();d=tri.GetOutput()
t=vtk_to_numpy(d.GetPolys().GetConnectivityArray()).reshape(-1,3);v=vtk_to_numpy(d.GetPoints().GetData())[t]
xyz=v.mean(axis=1);U=vtk_to_numpy(d.GetCellData().GetArray('U'));un=-U[:,2]
fig,axs=plt.subplots(2,3,figsize=(10,6.8));ports=[p for p in M['ports'] if p['region']=='Water_Xminus' and p['side']=='OUT']
for ax,port in zip(axs.ravel(),ports):
 center=np.array([port['x_mm'],port['y_mm']])*.001
 mask=np.linalg.norm(xyz[:,:2]-center,axis=1)<.00155
 points=(v[mask,:,:2]-center)*1000
 coll=PolyCollection(points,array=un[mask],cmap='RdBu_r',clim=(-4,4),edgecolors='none',rasterized=True);ax.add_collection(coll)
 ax.set(xlim=(-1.6,1.6),ylim=(-1.6,1.6),aspect='equal',title=port['patch'],xlabel='Local X [mm]',ylabel='Local Y [mm]')
fig.colorbar(coll,ax=axs.ravel().tolist(),fraction=.025,pad=.04,label='Velocity towards outlet (downwards) [m/s]')
fig.suptitle('Baseline iteration 500: 1 mm below outlet-pipe entrance (Z=-1 mm)\nBlue = locally upwards inside pipe; not reverse flow at the exit boundary',fontsize=11)
fig.subplots_adjust(top=.86,hspace=.42,wspace=.38,right=.88);fig.savefig(R/'diagnosis/internal-pipe-reverse-velocity.png',dpi=180);plt.close(fig)
