"""Geometry-only context of CAD OUT01..06 and divider movements."""
from pathlib import Path
import os,csv,math,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-seminar-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.font_manager import FontProperties
O=Path('/workspace/TBM/seminar-prep-20261125');CAD=O.parent/'cad-proposals-20261120'
ports=[p for p in csv.DictReader((CAD/'ports.csv').open()) if p['side']=='OUT' and p['region']=='Water_Xminus']
centre=-23.199977677315708;models=['B00','B01','C01','C02'];rows=[]
plt.rcParams['font.family']=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc').get_name()
fig,axs=plt.subplots(1,4,figsize=(13,6.2));yy=np.linspace(0,43,600)
for ax,n in zip(axs,models):
 offset=.5 if n=='B01' else 0;a=math.radians(-1 if n=='C01' else 1 if n=='C02' else 0)
 xx=centre+offset-np.tan(a)*yy;half=.5/math.cos(a)
 baseok=yy**2+(abs(centre)+.5)**2<47.9**2;ok=yy**2+np.maximum(abs(xx-half),abs(xx+half))**2<47.9**2
 ax.fill_betweenx(yy,centre-.5,centre+.5,where=baseok,color='#b4bcc3',alpha=.6)
 ax.plot(np.full_like(yy[baseok],centre-.5),yy[baseok],color='#6c7b88',lw=.8,ls='--')
 ax.fill_betweenx(yy,xx-half,xx+half,where=ok,color='#7d8994' if n=='B00' else '#e78a2d',alpha=.7)
 theta=np.linspace(np.pi/2,np.pi,500);ax.plot(47.9*np.cos(theta),47.9*np.sin(theta),c='#bac3ca',lw=.8)
 for p in ports:
  x,y=float(p['x_mm']),float(p['y_mm']);ax.add_patch(Circle((x,y),1.5,ec='#1676b9',fc='#e2f0f8',lw=1.1));ax.text(x-1,y+2.2,p['port'],fontsize=8,color='#1676b9',ha='center')
  delta=offset-math.tan(a)*y
  rows.append(dict(model=n,CAD_port=p['port'],pipe_x_mm=x,pipe_y_mm=y,divider_center_x_at_same_y_mm=centre+delta,divider_center_shift_from_B00_mm=delta,not_velocity_or_flow_result=True))
 ax.set(xlim=(-49,-20.3),ylim=(0,45),xlabel='X [mm]',title=n);ax.set_aspect('equal');ax.grid(alpha=.1)
axs[0].set_ylabel('Y [mm]：出口側')
fig.suptitle('左領域の6出口と隔壁：OUT_06付近で角度変更の移動量が大きい',fontsize=15,color='#14324b')
fig.text(.5,.015,'同じXY投影・縮尺。灰色/破線：B00、橙：変更後。水路の実断面や流量結果ではない。番号はCAD座標による。',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.055,1,.94));fig.savefig(O/'figures/outlet01-to06-geometry-context.png',dpi=210)
with (O/'geometry-only-outlet-divider-offsets.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Geometry-only outlet context and offsets written; no CFD')
