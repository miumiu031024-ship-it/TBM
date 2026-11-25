"""Export only flagged mesh cells, WITHOUT fields, and locate them geometrically."""
from pathlib import Path
import re,json,subprocess,shutil,math,csv,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-seminar-mpl')
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
R=Path(__file__).resolve().parent;O=Path('/workspace/TBM/seminar-prep-20261125');D=O/'mesh-quality';D.mkdir(exist_ok=True)
allrows=[];summaries=[]
for n in ['B00','C01','C02']:
 C=R/n;log=(C/'log.checkMesh-additional').read_text();m=re.search(r'Cells with small determinant .*?number of cells:\s*(\d+)',log);count=int(m.group(1)) if m else 0
 result=dict(model=n,small_determinant_cell_count=count,flow_solver_executed=False)
 if count:
  source=C/'VTK-quality/underdeterminedCells_0/internal.vtu'
  if not source.exists():
   cmd=f'source /usr/lib/openfoam/openfoam2412/etc/bashrc; cd /work/{n}; foamToVTK -constant -noZero -cellSet underdeterminedCells -no-fields -no-boundary -noFunctionObjects -with-ids -name VTK-quality'
   with (C/'log.foamToVTK-quality').open('w') as f:
    subprocess.run(['docker','exec','tbm-seminar-1125','/bin/bash','-c',cmd],stdout=f,stderr=subprocess.STDOUT,check=True)
  rr=vtk.vtkXMLUnstructuredGridReader();rr.SetFileName(str(source));rr.Update();data=rr.GetOutput();assert data.GetNumberOfCells()==count
  assert [data.GetCellData().GetArrayName(i) for i in range(data.GetCellData().GetNumberOfArrays())]==['cellID']
  cc=vtk.vtkCellCenters();cc.SetInputData(data);cc.Update();xyz=vtk_to_numpy(cc.GetOutput().GetPoints().GetData())*1000
  ids=vtk_to_numpy(data.GetCellData().GetArray('cellID'))
  cellset=C/'constant/polyMesh/sets/underdeterminedCells';raw=cellset.read_bytes();mm=re.search(rb'\n(\d+)\s*\(([^)]*)\)',raw)
  assert mm and int(mm.group(1))==count
  native_ids=[int(x) for x in mm.group(2).split()];assert len(native_ids)==count
  result['native_cell_ids_unordered']=native_ids
  result['subset_ID_warning']='VTK cellID is the exported subset index, NOT an asserted mapping to native global IDs. Native cellSet is retained separately.'
  rows=[]
  for i,p in zip(ids,xyz):
   plane_dist=[];normals=[]
   for c in [-23.199977677315708,23.199999069163503]:
    base=1 if n=='C01' else -1 if n=='C02' else 0;theta=math.radians(base if c>0 else -base)
    signed=(p[0]-c)*math.cos(theta)+p[1]*math.sin(theta);normals.append(signed);plane_dist.append(abs(abs(signed)-.5))
   region='Water_Xminus' if normals[0]<-.5 else 'Water_Xplus' if normals[1]>.5 else 'Water_Center'
   radius=float(np.linalg.norm(p));radial=min(abs(radius-44.9),abs(radius-47.9))
   row=dict(model=n,VTK_subset_index=int(i),x_mm=float(p[0]),y_mm=float(p[1]),z_mm=float(p[2]),radius_mm=radius,
    distance_to_nearest_sphere_mm=radial,distance_to_nearest_divider_face_mm=min(plane_dist),region_inferred_from_CAD_planes=region,
    center_definition='VTK cell representative point; not a solved physical quantity')
   rows.append(row);allrows.append(row)
  result['representative_points']=rows
  result['all_points_within_0p2mm_of_divider_and_sphere']=all(q['distance_to_nearest_sphere_mm']<.2 and q['distance_to_nearest_divider_face_mm']<.2 for q in rows)
  result['inferred_region_counts']={reg:sum(q['region_inferred_from_CAD_planes']==reg for q in rows) for reg in ['Water_Xminus','Water_Center','Water_Xplus']}
  shutil.copy2(source,D/f'{n}-small-determinant-cells.vtu');shutil.copy2(cellset,D/f'{n}-native-underdeterminedCells')
 summaries.append(result)
 print(n,'flagged cells',count,'root-proximity',result.get('all_points_within_0p2mm_of_divider_and_sphere'),flush=True)
(D/'small-determinant-cell-location-summary.json').write_text(json.dumps(summaries,indent=2))
if allrows:
 with (D/'small-determinant-cell-locations.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 plt.rcParams['font.family']=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc').get_name()
 fig=plt.figure(figsize=(11,5.6));phi=np.linspace(0,np.pi,200)
 for j,n in enumerate(['C01','C02']):
  ax=fig.add_subplot(1,2,j+1,projection='3d')
  for c in [-23.199977677315708,23.199999069163503]:
   base=1 if n=='C01' else -1;theta=math.radians(base if c>0 else -base);normal=np.array([math.cos(theta),math.sin(theta),0]);tangent=np.array([-math.sin(theta),math.cos(theta),0])
   for sign in [-1,1]:
    d=c*math.cos(theta)+sign*.5
    for radius in [44.9,47.9]:
     rc=math.sqrt(radius*radius-d*d);line=d*normal+rc*np.cos(phi)[:,None]*tangent+rc*np.sin(phi)[:,None]*np.array([0,0,1]);ax.plot(*line.T,c='#acb7c1',lw=.7)
  th=np.linspace(0,2*np.pi,300);ax.plot(47.9*np.cos(th),47.9*np.sin(th),np.zeros_like(th),c='#c9d0d6',lw=.8)
  rows=[r for r in allrows if r['model']==n]
  if rows:
   pp=np.array([[q['x_mm'],q['y_mm'],q['z_mm']] for q in rows]);ax.scatter(*pp.T,c='#d64b36',s=34,depthshade=False)
  ax.set(xlim=(-50,50),ylim=(-50,50),zlim=(0,51),xlabel='X [mm]',ylabel='Y [mm]',zlabel='Z [mm]',title=f'{n}：該当{len(rows)}セル');ax.view_init(22,-60);ax.set_box_aspect((100,100,51));ax.set_proj_type('ortho')
 fig.suptitle('追加検査で見つかった小さい行列式のセルの位置',fontsize=16,color='#14324b')
 fig.text(.5,.015,'赤点：メッシュセルの代表点。灰線：既存CAD寸法に基づく隔壁と球面の交線。流れの計算結果ではない。',ha='center',fontsize=9)
 fig.tight_layout(rect=(0,.07,1,.94));fig.savefig(O/'figures/quality-cell-locations.png',dpi=210)
