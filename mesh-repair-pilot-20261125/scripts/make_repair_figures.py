from pathlib import Path
import os,json,csv
os.environ['MPLCONFIGDIR']='/tmp/tbm-repair-mpl'
import numpy as np,vtk
from vtk.util.numpy_support import vtk_to_numpy,numpy_to_vtk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
R=Path(__file__).resolve().parent;O=Path('/workspace/TBM/mesh-repair-pilot-20261125');F=O/'figures';F.mkdir(exist_ok=True)
native_after={(q['model'],int(q['native_global_cell_id'])):float(q['checkMesh_internal_face_determinant']) for q in csv.DictReader((O/'repaired-target-cell-indicators.csv').open())} if (O/'repaired-target-cell-indicators.csv').exists() else {}
plt.rcParams['font.family']=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc').get_name()
fig=plt.figure(figsize=(12,5.5));axs=[]
for j,n in enumerate(['C01','C02']):
 changes=json.loads((R/(n+'-R07-point-adjustments.json')).read_text());ax=fig.add_subplot(1,2,j+1,projection='3d');axs.append(ax)
 data=vtk.vtkXMLUnstructuredGridReader();data.SetFileName('/workspace/TBM/seminar-prep-20261125/mesh-quality/'+n+'-small-determinant-cells.vtu');data.Update();grid=vtk.vtkUnstructuredGrid();grid.DeepCopy(data.GetOutput());xyz=vtk_to_numpy(grid.GetPoints().GetData());nativeids=np.full(6,-1,dtype=np.int64);nativeDet=np.zeros(6)
 for change in changes:
  old=np.array(change['original_position_m']);distance=np.linalg.norm(xyz-old,axis=1);i=int(np.argmin(distance));assert distance[i]<1e-8
  cells=vtk.vtkIdList();grid.GetPointCells(i,cells);assert cells.GetNumberOfIds()==1;cell=int(cells.GetId(0));nativeids[cell]=change['original_bad_cell'];nativeDet[cell]=native_after.get((n,change['original_bad_cell']),change['local_after_min_determinant']);grid.GetPoints().SetPoint(i,change['new_position_m'])
 assert (nativeids>=0).all();a=numpy_to_vtk(nativeids);a.SetName('nativeGlobalCellID');grid.GetCellData().AddArray(a);a=numpy_to_vtk(nativeDet);a.SetName('postRepairDeterminant');grid.GetCellData().AddArray(a);grid.GetPoints().Modified()
 w=vtk.vtkXMLUnstructuredGridWriter();w.SetFileName(str(O/(n+'-six-repaired-target-cells.vtu')));w.SetInputData(grid);w.Write()
 p=np.array([q['original_position_m'] for q in changes])*1000;ax.scatter(*p.T,c='#bf4433',s=30)
 for p0,q in zip(p,changes):ax.text(*p0,str(q['original_bad_cell']),fontsize=6)
 for x in [-23.2,23.2]:
  th=np.linspace(0,np.pi,100);yaw=np.radians((1 if n=='C01' else -1)*(1 if x>0 else -1));normal=np.array([np.cos(yaw),np.sin(yaw),0]);tangent=np.array([-np.sin(yaw),np.cos(yaw),0]);d=x*np.cos(yaw);rad=np.sqrt(47.9**2-d*d);line=d*normal+rad*np.cos(th)[:,None]*tangent+rad*np.sin(th)[:,None]*np.array([0,0,1]);ax.plot(*line.T,c='#a8b5c0',lw=.8)
 th=np.linspace(0,2*np.pi,100);ax.plot(47.9*np.cos(th),47.9*np.sin(th),np.zeros(100),c='#a8b5c0',lw=.8)
 ax.set(xlim=(-50,50),ylim=(-50,50),zlim=(-2,50),xlabel='X [mm]',ylabel='Y [mm]',zlabel='Z [mm]',title=n+'：内部頂点6点のみ修正');ax.view_init(22,-60);ax.set_box_aspect((100,100,52));ax.set_proj_type('ortho')
fig.suptitle('小さい行列式が生じたセルの位置と修正対象',fontsize=16,color='#14324b');fig.text(.5,.015,'赤点は移動した内部頂点の元の位置。数µmの移動は実寸表示では識別困難。灰線はCAD寸法に基づく概略輪郭。境界メッシュ・CADは変更していない。',ha='center',fontsize=10);fig.tight_layout(rect=(0,.05,1,.93));fig.savefig(F/'repair-locations.png',dpi=180);plt.close(fig)
# Determinant and displacement plot, no solved-flow images.
fig,ax=plt.subplots(1,2,figsize=(11,4.5));raw=list(csv.DictReader((O/'original-bad-cell-indicators.csv').open()));labels=[];before=[];after=[];distance=[]
for n in ['C01','C02']:
 for q in json.loads((R/(n+'-R07-point-adjustments.json')).read_text()):
  labels.append(n+'\n'+str(q['original_bad_cell']));before.append(float(next(r['checkMesh_internal_face_determinant'] for r in raw if r['model']==n and int(r['native_global_cell_id'])==q['original_bad_cell'])));after.append(native_after.get((n,q['original_bad_cell']),q['local_after_min_determinant']));distance.append(np.linalg.norm(q['displacement_um']))
x=np.arange(12);ax[0].scatter(x,before,label='修正前',color='#ba493b');ax[0].scatter(x,after,label='修正後（native面から再計算）',color='#247a76');ax[0].axhline(.001,color='#333',ls='--',lw=1,label='検査基準 0.001');ax[0].set_yscale('log');ax[0].set_ylabel('checkMeshの行列式指標');ax[0].legend(fontsize=8);ax[1].bar(x,distance,color='#247a76');ax[1].set_ylabel('内部頂点の移動距離 [µm]')
for a in ax:a.set_xticks(x,labels,fontsize=6,rotation=60);a.grid(axis='y',alpha=.2)
fig.suptitle('修正は境界形状を保ったメッシュ内部の調整',fontsize=15);fig.tight_layout();fig.savefig(F/'determinant-and-displacement.png',dpi=180)
