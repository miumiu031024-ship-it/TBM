"""Scientific CAD figures; all geometry from existing STEP. No CFD colors."""
from pathlib import Path
import os,json,csv,math
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-seminar-mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-seminar-cache')
import numpy as np
import cadquery as cq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch,Circle,Rectangle,FancyArrowPatch,FancyBboxPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from OCP.BRepTools import BRepTools
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRep import BRep_Tool
from OCP.TopLoc import TopLoc_Location

ROOT=Path(__file__).resolve().parent
OUT=Path('/workspace/TBM/seminar-prep-20261125');(OUT/'figures').mkdir(exist_ok=True)
O=OUT/'figures';CAD=OUT.parent/'cad-proposals-20261120'
font='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
plt.rcParams.update({'font.family':FontProperties(fname=font).get_name(),'font.size':10,
 'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
models={'B00':'B00_baseline','B01':'B01_partition_inward_050','C01':'C01_partition_yaw_plus1','C02':'C02_partition_yaw_minus1'}
water={n:cq.importers.importStep(str(CAD/'water'/f'{f}_water.step')).val() for n,f in models.items()}
metal={n:cq.importers.importStep(str(CAD/'partition-reference'/f'{f}_partition_reference.step')).val() for n,f in models.items()}
ports=list(csv.DictReader((CAD/'ports.csv').open()))
colors=['#3c9fba','#75b886','#9a86be'];orange='#e78a2d';ink='#14324b'
def parts(s):return sorted(s.Solids(),key=lambda a:a.Center().x)
def mesh(s,tol=.05):
 triangles=[]
 for solid in s.Solids() or [s]:
  BRepTools.Clean_s(solid.wrapped);BRepMesh_IncrementalMesh(solid.wrapped,tol,False,.18)
  for f in solid.Faces():
   loc=TopLoc_Location();p=BRep_Tool.Triangulation_s(f.wrapped,loc)
   if p is None:continue
   v=np.array([p.Node(i).Transformed(loc.Transformation()).Coord() for i in range(1,p.NbNodes()+1)])
   ind=np.array([p.Triangle(i).Get() for i in range(1,p.NbTriangles()+1)])-1
   triangles.extend(v[ind])
 return np.asarray(triangles)
cache={n:[mesh(s) for s in parts(water[n])] for n in models}
mcache={n:mesh(s) for n,s in metal.items()}
def setup3d(ax):
 ax.set(xlim=(-51,51),ylim=(-51,51),zlim=(-12,52));ax.set_box_aspect((102,102,64));ax.view_init(22,-60);ax.set_axis_off();ax.set_proj_type('ortho')
def draw3d(ax,name,domain_colors=True):
 for t,c in zip(cache[name],colors):
  ax.add_collection3d(Poly3DCollection(t,facecolor=c if domain_colors else '#80b5c6',alpha=.28,edgecolor='none',rasterized=True))
 ax.add_collection3d(Poly3DCollection(mcache[name],facecolor='#8a939c' if name=='B00' else orange,alpha=.92,edgecolor='none',rasterized=True))
 setup3d(ax)
def save(fig,name):
 fig.savefig(O/name,dpi=210,bbox_inches='tight');plt.close(fig);print(name,flush=True)

fig=plt.figure(figsize=(13,6));ax=fig.add_subplot(121,projection='3d');draw3d(ax,'B00')
ax.set_title('B00：球面水路＋36本の短管',color=ink,pad=0,fontsize=15)
ax.text2D(.03,.03,'半径 44.9–47.9 mm（水路厚3 mm）\n管内径3 mm／直管長10 mm',transform=ax.transAxes,fontsize=10)
ax=fig.add_subplot(122)
th=np.linspace(0,2*np.pi,500);ax.plot(47.9*np.cos(th),47.9*np.sin(th),c='#bac3ca',lw=1)
for x in [-23.1999776773,23.1999990692]:
 ymax=math.sqrt(47.9**2-x*x);ax.add_patch(Rectangle((x-.5,-ymax),1,2*ymax,color='#9099a2',alpha=.6))
for side,c in [('IN','#1676b9'),('OUT','#d64b49')]:
 for p in [p for p in ports if p['side']==side]:
  x,y=float(p['x_mm']),float(p['y_mm']);ax.add_patch(Circle((x,y),1.5,fc=c,ec='white',lw=.4))
  factor=1.22 if abs(x)<15 and int(p['port'][-2:])%2==0 else 1.17
  ax.text(x*factor,y*factor,p['port'],ha='center',va='center',fontsize=6.4,color=c)
ax.text(-37,0,'左側\n6入口＋6出口',ha='center',va='center',color=colors[0]);ax.text(0,0,'中央\n6入口＋6出口',ha='center',va='center',color=colors[1]);ax.text(37,0,'右側\n6入口＋6出口',ha='center',va='center',color=colors[2])
ax.set(xlim=(-60,60),ylim=(-60,60),xlabel='X [mm]',ylabel='Y [mm]',title='36開口のXY配置：入口Y<0／出口Y>0');ax.set_aspect('equal');ax.grid(alpha=.1)
fig.suptitle('既存STEPの確認：3つの独立した水領域・隔壁2枚',fontsize=18,color=ink,y=1.01)
fig.text(.5,.01,'形状図：色は水領域の区別。速度分布ではない。灰色は水路内の隔壁参考部で、金属全体のCADではない。',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.05,1,.96));save(fig,'B00-overview-and-ports.png')

def section(ax,shape,plane,color,axes):
 cut=shape.intersect(plane)
 for f in cut.Faces():
  vv,tt=f.tessellate(.005,.12)
  if tt:
   v=np.array([p.toTuple() for p in vv]);t=v[np.array(tt,dtype=int)]
   ax.add_collection(PolyCollection(t[:,:,axes],fc=color,ec='none'))
 for e in cut.Edges():
  n=max(2,min(250,int(e.Length()/.08)+1));v=np.array([e.positionAt(t).toTuple() for t in np.linspace(0,1,n)])
  ax.plot(v[:,axes[0]],v[:,axes[1]],c=ink,lw=.45)
fig,axs=plt.subplots(1,2,figsize=(13,5.5),gridspec_kw={'width_ratios':[1.7,1]})
plane=cq.Face.makePlane(125,125,cq.Vector(0,0,22),cq.Vector(0,1,0))
for s,c in zip(parts(water['B00']),colors):
 section(axs[0],s,plane,c,[0,2]);section(axs[1],s,plane,c,[0,2])
for ax in axs:section(ax,metal['B00'],plane,'#929aa3',[0,2]);ax.set_aspect('equal');ax.set_xlabel('X [mm]');ax.set_ylabel('Z [mm]');ax.grid(alpha=.1)
axs[0].set(xlim=(-51,51),ylim=(-2,51),title='Y=0の実CAD断面（この面は管を横切らない）')
axs[0].text(0,21,'内側の非水領域\n模型の内壁・支持構造は別途設計',ha='center',va='center',fontsize=12,color='#6d7a87')
axs[0].annotate('水路の径方向厚さ3 mm',xy=(-34,32),xytext=(-48,15),arrowprops={'arrowstyle':'->','color':ink},fontsize=10)
axs[1].set(xlim=(20.5,26),ylim=(36.3,44),title='右隔壁付け根の拡大：B00は鋭い角')
axs[1].annotate('隔壁相当幅1 mm',xy=(23.2,40),xytext=(20.7,43),arrowprops={'arrowstyle':'->','color':ink},fontsize=10)
fig.suptitle('B00の断面：水の空間と、製作用の壁・構造体を区別する',color=ink,fontsize=17)
fig.text(.5,.01,'既存STEPは水領域。透明外壁の厚さ、締結部、シール、継手は含まれていない。',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.05,1,.94));save(fig,'B00-real-CAD-section.png')

labels={'B00':'基準：隔壁位置・角度を維持','B01':'中央へ各0.50 mm平行移動','C01':'出口側で左右を広げる方向','C02':'出口側で左右を狭める方向'}
fig=plt.figure(figsize=(16,4.5))
for i,n in enumerate(models):
 ax=fig.add_subplot(1,4,i+1,projection='3d');draw3d(ax,n,False)
 ax.set_title(n+'\n'+labels[n],fontsize=12,color=ink,pad=0)
fig.suptitle('同じ視点・縮尺での形状比較',fontsize=20,color=ink,y=1.03)
fig.text(.5,.02,'水領域：水色／基準隔壁：灰色／変更した隔壁：橙色。球面外形・36管・3領域を維持。CFD結果ではない。',ha='center',fontsize=11)
fig.subplots_adjust(left=0,right=1,top=.88,bottom=.06,wspace=0);save(fig,'four-models-same-view.png')

fig,axs=plt.subplots(1,4,figsize=(15,5))
base_c=-23.199977677315708;yy=np.linspace(33,43,400);pipe=next(p for p in ports if p['port']=='OUT_06');px,py=float(pipe['x_mm']),float(pipe['y_mm'])
for ax,n in zip(axs,models):
 shift=.5 if n=='B01' else 0;theta=np.deg2rad(-1 if n=='C01' else 1 if n=='C02' else 0)
 xx=base_c+shift-np.tan(theta)*yy;half=.5/np.cos(theta)
 base_ok=yy**2+max(abs(base_c-.5),abs(base_c+.5))**2<=47.9**2
 current_ok=yy**2+np.maximum(abs(xx-half),abs(xx+half))**2<=47.9**2
 ax.fill_betweenx(yy,base_c-.5,base_c+.5,where=base_ok,color='#acb4bc',alpha=.5)
 ax.plot(np.full_like(yy[base_ok],base_c-.5),yy[base_ok],c='#6e7883',ls='--',lw=.9);ax.plot(np.full_like(yy[base_ok],base_c+.5),yy[base_ok],c='#6e7883',ls='--',lw=.9)
 if n!='B00':ax.fill_betweenx(yy,xx-half,xx+half,where=current_ok,color=orange,alpha=.72)
 else:ax.fill_betweenx(yy,xx-half,xx+half,where=current_ok,color='#7c8995',alpha=.65)
 ax.add_patch(Circle((px,py),1.5,ec='#1676b9',fc='#dfedf8',lw=1.7));ax.plot(px,py,'+',c='#1676b9')
 ax.text(px,py-2.15,'OUT_06\n位置・内径は同じ',ha='center',fontsize=9,color='#1676b9')
 if n!='B00':
  delta=shift-np.tan(theta)*py
  ax.annotate('',xy=(base_c+delta,py),xytext=(base_c,py),arrowprops={'arrowstyle':'<->','color':'#ac560a','lw':1.2})
  ax.text(-24.6,41.65,f'このYでの中心移動\n{delta:+.2f} mm',ha='center',fontsize=9,color='#ac560a')
 ax.set(xlim=(-29,-20.5),ylim=(33,43),xlabel='X [mm]',title=n+'\n'+labels[n]);ax.set_aspect('equal');ax.grid(alpha=.13)
axs[0].set_ylabel('Y [mm]：上側が出口側')
fig.suptitle('出口付近の隔壁位置の差：OUT_06付近のXY投影・等縮尺',fontsize=18,color=ink,y=1.02)
fig.text(.5,.01,'破線・灰色：B00の隔壁。橙色：変更後。球面水路と短管の投影図であり、水平面の流体断面ではない。番号はCAD座標による。',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.12,1,.93));save(fig,'divider-outlet-detail.png')

fig,axs=plt.subplots(2,1,figsize=(13,5.8));
for ax in axs:ax.set(xlim=(0,13),ylim=(0,2.6));ax.axis('off')
def box(ax,x,y,w,h,label,color):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.04,rounding_size=.08',fc=color,ec='none'))
 ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=11,color=ink)
def arrow(ax,x0,x1,y):ax.annotate('',xy=(x1,y),xytext=(x0,y),arrowprops={'arrowstyle':'->','lw':1.7,'color':ink})
ax=axs[0];ax.text(0,2.35,'現模型：ユーザーの説明に基づく流れの接続関係（寸法・位置は未確認）',fontsize=13,color=ink)
box(ax,.1,.8,2.5,.95,'追加の接続管\n15本','#dcebf6');box(ax,3.25,.8,4,.95,'大きな共通空間\n空気残留・気泡の上部集積を観察','#ffe5bf');box(ax,8,.8,4.5,.95,'冷却流路\n共通空間との接続詳細は未確認','#e3edf0');arrow(ax,2.7,3.15,1.28);arrow(ax,7.35,7.9,1.28)
ax.text(.1,.25,'設置：球面が上。下方の管から水・気泡を注入。流路厚・隔壁・18入口/18出口との対応・排気口は要確認。',fontsize=9.5,color='#6d7a87')
ax=axs[1];ax.text(0,2.35,'STEPに近づける案：観察する水路と、給排水・排気の設計を分ける',fontsize=13,color=ink)
box(ax,.1,.8,2.5,.95,'給水・分配系\n実測流量を確認','#dcebf6');box(ax,3.25,.8,5.5,.95,'18入口 → 3独立球面水路 → 18出口\n追加の共通空間は観察部に設けない','#d8eee5');box(ax,9.4,.8,3.1,.95,'排水・集合系\n境界への影響を確認','#e3edf0');arrow(ax,2.7,3.15,1.28);arrow(ax,8.85,9.3,1.28)
ax.text(.1,.25,'必要：透明内外壁／外周継手／シール／真の高所での排気。可視化用気泡の注入は排気と分ける。',fontsize=10,color='#6d7a87')
fig.suptitle('透明模型の見直し：追加空間の除去だけでは、気泡の浮上は解消しない',fontsize=18,color=ink,y=1.01)
fig.tight_layout(rect=(0,0,1,.95));save(fig,'experiment-model-concept.png')
(OUT/'figure-provenance.json').write_text(json.dumps(dict(source='Existing STEP files only; no geometry altered',models=models,camera=dict(elevation_deg=22,azimuth_deg=-60,projection='orthographic',xlim_mm=[-51,51],ylim_mm=[-51,51],zlim_mm=[-12,52]),figures_not_CFD=True,current_experiment_schematic_source='User description: dome up, bottom tube water/bubble inlet; dimensions and detailed connections not assumed'),indent=2,ensure_ascii=False))
