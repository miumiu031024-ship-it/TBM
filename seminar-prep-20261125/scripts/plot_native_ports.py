from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-seminar-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle
from matplotlib.font_manager import FontProperties
O=Path('/workspace/TBM/seminar-prep-20261125')
plt.rcParams['font.family']=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc').get_name()
fig,axs=plt.subplots(1,3,figsize=(12,4.7))
for ax,n in zip(axs,['B00','C01','C02']):
 q=json.loads((O/f'mesh-quality/{n}-OUT06-native-mesh.json').read_text())
 ax.add_collection(PolyCollection(q['polygons_local_mm'],fc='#f0f5f8',ec='#314e63',lw=.5))
 ax.add_patch(Circle((0,0),1.5,fill=False,ec='#d67a25',ls='--',lw=1.4))
 ax.set(xlim=(-1.7,1.7),ylim=(-1.7,1.7),xlabel='管中心からのX [mm]',ylabel='管中心からのY [mm]',title=f'{n}：OUT_06端面\n{q["native_faces"]}面');ax.set_aspect('equal');ax.grid(alpha=.1)
fig.suptitle('実際に生成した出口端面のメッシュ：同じ位置・縮尺',fontsize=16,color='#14324b')
fig.text(.5,.01,'橙破線：CAD内径3 mm。幾何・格子の確認図であり、速度や圧力の計算結果ではない。',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.14,1,.91));fig.savefig(O/'figures/native-outlet-mesh-check.png',dpi=210)
