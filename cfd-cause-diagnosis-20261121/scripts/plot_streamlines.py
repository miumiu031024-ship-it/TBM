"""Streamlines in one actual saved field; no physical-time interpretation."""
from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkFiltersFlowPaths import vtkStreamTracer
from vtkmodules.util.numpy_support import vtk_to_numpy,numpy_to_vtk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
R=Path(__file__).resolve().parent;M=json.loads((R/'SeedSide/case-metadata.json').read_text())
r=vtkXMLUnstructuredGridReader();r.SetFileName(str(R/'SeedSide/VTK/SeedSide_0/internal.vtu'));r.Update();d=r.GetOutput();assert d.GetPointData().GetArray('U') is not None
d.GetPointData().SetActiveVectors('U')
ports=[p for p in M['ports'] if p['region']=='Water_Xminus' and p['side']=='IN'];seeds=[]
for p in ports:
 for i in range(12):
  rr=.001*np.sqrt((i+.5)/12);ang=i*np.pi*(3-np.sqrt(5));seeds.append((p['x_mm']*.001+rr*np.cos(ang),p['y_mm']*.001+rr*np.sin(ang),-.0085))
pts=vtkPoints();pts.SetData(numpy_to_vtk(np.array(seeds),deep=True));seed=vtkPolyData();seed.SetPoints(pts)
s=vtkStreamTracer();s.SetInputData(d);s.SetSourceData(seed);s.SetIntegratorTypeToRungeKutta45();s.SetIntegrationDirectionToForward();s.SetIntegrationStepUnit(s.LENGTH_UNIT);s.SetInitialIntegrationStep(.0001);s.SetMinimumIntegrationStep(.00001);s.SetMaximumIntegrationStep(.0002);s.SetMaximumPropagation(.30);s.SetMaximumNumberOfSteps(10000);s.SetTerminalSpeed(1e-5);s.Update();out=s.GetOutput()
p=vtk_to_numpy(out.GetPoints().GetData())*1000;u=vtk_to_numpy(out.GetPointData().GetArray('U'));speed=np.linalg.norm(u,axis=1)
segments=[];colors=[];ends=[]
for i in range(out.GetNumberOfLines()):
 line=out.GetCell(i);ids=np.array([line.GetPointId(j) for j in range(line.GetNumberOfPoints())]);q=p[ids]
 segments.extend(np.stack([q[:-1],q[1:]],axis=1));colors.extend((speed[ids][:-1]+speed[ids][1:])/2);ends.append(q[-1])
fig=plt.figure(figsize=(10,6));ax=fig.add_subplot(111,projection='3d')
lc=Line3DCollection(segments,array=np.array(colors),cmap='turbo',norm=plt.Normalize(0,4),linewidth=.6);ax.add_collection(lc)
for p0 in [p for p in M['ports'] if p['region']=='Water_Xminus']:
 ax.plot([p0['x_mm']]*2,[p0['y_mm']]*2,[-10,0],c='#aaa',lw=1)
x=-23.6999776773;rr=np.sqrt(46.4**2-x*x);a=np.linspace(0,np.pi,120);ax.plot(x+np.zeros_like(a),rr*np.cos(a),rr*np.sin(a),c='#c55',lw=1,label='Divider at X=-23.7 mm')
ax.set(xlim=(-49,-21),ylim=(-43,43),zlim=(-10,43),xlabel='X [mm]',ylabel='Y [mm]',zlabel='Z [mm]');ax.set_box_aspect((28,86,53));ax.view_init(22,-145);ax.legend(loc='upper left',fontsize=7)
fig.colorbar(lc,ax=ax,shrink=.55,label='Speed [m/s]');fig.suptitle('Baseline left domain: streamlines from inlet-pipe seeds, iteration 500\nUnsettled saved field; not particle tracks, not evidence of physical periodicity',fontsize=11)
fig.savefig(R/'diagnosis/baseline-streamlines.png',dpi=180,bbox_inches='tight');plt.close(fig)
(R/'diagnosis/streamline-metadata.json').write_text(json.dumps(dict(seed_count=len(seeds),integrated_lines=out.GetNumberOfLines(),maximum_propagation_m=.30,actual_points=out.GetNumberOfPoints(),interpretation='Instantaneous streamlines of a saved unsettled steady field only'),indent=2))
print('Streamlines:',out.GetNumberOfLines(),'points:',out.GetNumberOfPoints())
