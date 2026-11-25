from pathlib import Path
import os,json,math
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-seminar-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
O=Path('/workspace/TBM/seminar-prep-20261125')
q=json.loads((O.parent/'cad-proposals-20261120/geometry-validation.json').read_text())
lo1,hi1,lo2,hi2=q['original_partition_face_x_mm']
plt.rcParams['font.family']=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc').get_name()
x=np.linspace(-47.9,47.9,5001);zo=np.sqrt(np.maximum(47.9**2-x*x,0));zi=np.sqrt(np.maximum(44.9**2-x*x,0))
fig,ax=plt.subplots(figsize=(11,5.4))
for mask,color in [(x<=lo1,'#3c9fba'),((x>=hi1)&(x<=lo2),'#75b886'),(x>=hi2,'#9a86be'),((x>lo1)&(x<hi1),'#9bA5ae'),((x>lo2)&(x<hi2),'#9bA5ae')]:
 ax.fill_between(x,zi,zo,where=mask,color=color)
ax.plot(x,zo,color='#486276',lw=.7);ax.plot(x,zi,color='#486276',lw=.7)
locations=[(lo1,math.sqrt(47.9**2-lo1**2),'左側の高所\nZ≈41.6 mm',(-43,49)),(0,47.9,'中央の高所\nZ=47.9 mm',(0,56)),(hi2,math.sqrt(47.9**2-hi2**2),'右側の高所\nZ≈41.6 mm',(43,49))]
for xx,zz,label,pos in locations:
 ax.plot(xx,zz,'o',ms=6,color='#c55d18');ax.annotate(label,xy=(xx,zz),xytext=pos,ha='center',fontsize=11,color='#a85312',arrowprops={'arrowstyle':'->','color':'#a85312','lw':1.3})
ax.text(0,22,'3領域は独立しているため、\n中央頂部の排気だけで左右も抜けるとは限らない。',ha='center',fontsize=13,color='#14324b')
ax.set(xlim=(-54,54),ylim=(-1,61),xlabel='X [mm]',ylabel='Z [mm]');ax.set_aspect('equal');ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.1)
ax.set_title('球面が上の設置：B00の各水領域の高所と排気の検討位置',fontsize=16,color='#14324b')
fig.text(.5,.015,'B00のY=0断面（水側寸法）に概念を追記。排気穴・継手は未作成。実験模型の排気位置を再現した図ではない。',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.055,1,1));fig.savefig(O/'figures/B00-high-point-vent-concept.png',dpi=210)
print('CAD-based high-point diagram written; no CAD altered')
