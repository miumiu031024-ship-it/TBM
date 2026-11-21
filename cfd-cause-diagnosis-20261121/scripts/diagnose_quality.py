from pathlib import Path
import re,json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkFiltersCore import vtkCellCenters
from vtkmodules.util.numpy_support import vtk_to_numpy
R=Path(__file__).resolve().parent;C=R/'SeedSide'
r=vtkXMLUnstructuredGridReader();r.SetFileName(str(C/'VTK/SeedSide_0/internal.vtu'));r.Update();d=r.GetOutput()
cc=vtkCellCenters();cc.SetInputData(d);cc.CopyArraysOff();cc.Update();xyz=vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
b=(C/'constant/polyMesh/sets/concaveCells').read_text();tail=b.split('// * * *')[1];mt=re.search(r'(\d+)\s*\(\s*([\d\s]+)\)',tail,re.S);assert mt
ids=np.fromstring(mt.group(2),dtype=np.int64,sep=' ');assert len(ids)==int(mt.group(1)) and len(ids)==14561 and ids.max()<len(xyz)
flag=np.zeros(len(xyz),bool);flag[ids]=True;radii=np.linalg.norm(xyz,axis=1)
regions={'All':np.ones(len(xyz),bool),'Straight_pipes_z_negative':xyz[:,2]<0,'Shell_near_pipe_junction_0_to_3mm':(xyz[:,2]>=0)&(xyz[:,2]<.003),'Near_divider_within_0p6mm':(abs(xyz[:,0]+.0236999776773)<.0006)&(xyz[:,2]>=0),'Shell_inner_0p3mm':(radii<.0452)&(xyz[:,2]>=0),'Shell_outer_0p3mm':(radii>.0476)&(xyz[:,2]>=0),'Shell_mid_gap_45p2_to_47p6mm':(radii>=.0452)&(radii<=.0476)&(xyz[:,2]>=0)}
a=[]
for n,m in regions.items():a.append(dict(region=n,cells=int(m.sum()),concave_cells=int((flag&m).sum()),concave_fraction=float(flag[m].mean())))
(R/'diagnosis/mesh-quality-location.json').write_text(json.dumps(dict(actual_cells=len(xyz),regions_overlap=True,regions=a,interpretation='Location audit only; correlation or concavity alone does not identify instability cause'),indent=2))
print(json.dumps(a,indent=2))
