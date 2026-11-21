"""Actual CAD section for explaining A01; no CFD contours in this figure."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-a01-plots')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-a01-cache')
import cadquery as cq
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch,Rectangle,ConnectionPatch
R=Path(__file__).resolve().parent;C=Path('/workspace/TBM/cad-proposals-20261120')
base=cq.importers.importStep(str(C/'water/B00_baseline_water.step')).val()
candidate=cq.importers.importStep(str(C/'water/A01_root_R025_water.step')).val()
metal=cq.importers.importStep(str(C/'partition-reference/B00_baseline_partition_reference.step')).val()
plane=cq.Face.makePlane(120,120,cq.Vector(0,0,20),cq.Vector(0,1,0))
sections=[s.intersect(plane) for s in [base,metal,base.cut(candidate)]]
water=candidate.intersect(plane)
def draw(ax):
    for s,col in zip(sections,['#b9deef','#999999','#ed8b20']):
        v,t=s.tessellate(.002,.04);v=np.array([p.toTuple() for p in v]);tri=v[np.array(t,dtype=int)]
        ax.add_collection(PolyCollection(tri[:,:,[0,2]],facecolor=col,edgecolor='none'))
    for e in water.Edges():
        n=max(2,min(3000,int(e.Length()/.02)+1));v=np.array([e.positionAt(t).toTuple() for t in np.linspace(0,1,n)])
        ax.plot(v[:,0],v[:,2],c='#1b4a63',lw=.7)
    ax.set_aspect('equal');ax.set(xlabel='X [mm]',ylabel='Z [mm]');ax.grid(alpha=.12)
fig,(ax,zoom)=plt.subplots(1,2,figsize=(11,6),gridspec_kw={'width_ratios':[1.5,1]})
draw(ax);draw(zoom)
ax.set(xlim=(-49,49),ylim=(-1,49),title='Location in Y=0 section')
ax.add_patch(Rectangle((21.2,36.9),4,6.1,fill=False,edgecolor='#bd3029',lw=1.2))
ax.annotate('Divider junction',xy=(23.2,40),xytext=(-5,26),arrowprops=dict(arrowstyle='->',color='#bd3029'),color='#7d2925')
zoom.set(xlim=(21.2,25.2),ylim=(36.9,43),title='A01: enlarged junction')
zoom.annotate('R0.25 mm',xy=(23.78,41.50),xytext=(23.1,42.55),arrowprops=dict(arrowstyle='->',color='#874c12'),color='#874c12',fontsize=11)
fig.legend(handles=[Patch(color='#b9deef',label='Water'),Patch(color='#999999',label='Original divider reference'),Patch(color='#ed8b20',label='Added metal / removed water')],loc='lower center',ncol=3,fontsize=10)
fig.suptitle('A01: small rounds at the divider / spherical-wall junctions\nR0.25 mm; water volume reduced by about 0.055%. CAD GEOMETRY ONLY',fontsize=13)
fig.tight_layout(rect=(0,.07,1,.89));fig.savefig(R/'A01-change-explained.png',dpi=180);plt.close(fig)
print('Rendered actual CAD section explaining A01.')
