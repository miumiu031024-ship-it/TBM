"""Localize measured changes, using actual B00 probe histories at 400/500."""
from pathlib import Path
import os,json,csv,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;O=R/'diagnosis';source=Path('/workspace/scratch/tbm-A01-cfd-20261120/results/B00/internal-iteration-monitor.npz')
d=np.load(source);ts=d['iterations'];ia=np.flatnonzero(ts==400)[0];ib=np.flatnonzero(ts==500)[0]
good=d['valid'][ia]&d['valid'][ib];xyz=d['points_m'][good];a=d['velocity_m_s'][ia,good];b=d['velocity_m_s'][ib,good];delta=np.linalg.norm(b-a,axis=1)
x=xyz[:,0];dist=np.min(abs(x[:,None]-np.array([-.0236999776773,-.0226999776773,.0226999990692,.0236999990692])[None,:]),axis=1)
regions={'All':np.ones(len(xyz),bool),'Within_1mm_of_divider_face':dist<.001,'Away_from_divider_faces_1mm':dist>=.001,'Lower_shell_z_under_5mm':xyz[:,2]<.005,'Mid_shell_z_5_to_20mm':(xyz[:,2]>=.005)&(xyz[:,2]<=.020),'Upper_shell_z_over_20mm':xyz[:,2]>.020}
rows=[]
for n,m in regions.items():
 rows.append(dict(region=n,points=int(m.sum()),mean_vector_change_m_s=float(delta[m].mean()),vector_relative_L2_change=float(np.sqrt(((b[m]-a[m])**2).sum()/(a[m]**2).sum())),fraction_points_changing_over_0p2_m_s=float((delta[m]>.2).mean())))
report=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),from_iteration=400,to_iteration=500,regions_overlap=True,regions=rows,interpretation='Steady numerical iteration difference; neither physical frequency nor causal attribution')
(O/'iteration-change-location.json').write_text(json.dumps(report,indent=2))
np.savez_compressed(O/'iteration-change-samples.npz',points_m=xyz,velocity_iteration400=a,velocity_iteration500=b,vector_difference_m_s=delta)
fig=plt.figure(figsize=(13,4.5))
for i,(val,title,maximum) in enumerate([(np.linalg.norm(a,axis=1),'Iteration 400 speed [m/s]',1.5),(np.linalg.norm(b,axis=1),'Iteration 500 speed [m/s]',1.5),(delta,'Vector difference magnitude [m/s]',1)]):
 ax=fig.add_subplot(1,3,i+1,projection='3d');p=xyz*1000
 sc=ax.scatter(p[:,0],p[:,1],p[:,2],c=val,s=3,cmap='turbo',vmin=0,vmax=maximum,depthshade=False,linewidth=0)
 ax.set(xlim=(-49,49),ylim=(-49,49),zlim=(0,49),title=title);ax.set_box_aspect((98,98,49));ax.view_init(22,-55);ax.set_axis_off();fig.colorbar(sc,ax=ax,shrink=.45,pad=0)
fig.suptitle('Baseline: actual common internal probes at r=46.4 mm\nSteady iteration drift, not a physical-time sequence',fontsize=11);fig.savefig(O/'where-internal-field-changes.png',dpi=180,bbox_inches='tight');plt.close(fig)
print(json.dumps(report,indent=2))
