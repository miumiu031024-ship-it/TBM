"""CAD geometry figures only. Colors are materials, never CFD velocity."""
import os,json,argparse,math
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-cad-cache')
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cad-mpl')
from pathlib import Path
import numpy as np
import cadquery as cq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
P=argparse.ArgumentParser();P.add_argument('folder',type=Path);args=P.parse_args();O=args.folder
D=json.loads((O/'geometry-validation.json').read_text())
models=D['models'];water={m['id']:cq.importers.importStep(str(O/m['water_file'])).val() for m in models};metal={m['id']:cq.importers.importStep(str(O/m['partition_reference_file'])).val() for m in models}
def mesh(s,tol=.08):
 vv,tt=s.tessellate(tol,.12);v=np.array([p.toTuple() for p in vv]);return v[np.array(tt,dtype=int)]
titles=['B00: unchanged baseline','A01: divider-root R0.25 mm','A02: divider-root R0.50 mm','B01: dividers inward0.50 mm','C01: mirrored yaw +1 deg','C02: mirrored yaw -1 deg']
fig=plt.figure(figsize=(13,8))
for j,m in enumerate(models):
 ax=fig.add_subplot(2,3,j+1,projection='3d')
 ax.add_collection3d(Poly3DCollection(mesh(water[m['id']]),facecolor='#79b7cf',alpha=.19,edgecolor='none',zsort='average'))
 ax.add_collection3d(Poly3DCollection(mesh(metal[m['id']],.03),facecolor='#e98524',edgecolor='none',alpha=1,zsort='average'))
 ax.set(xlim=(-51,51),ylim=(-51,51),zlim=(-12,52));ax.set_box_aspect((102,102,64));ax.view_init(22,-55);ax.set_axis_off();ax.set_title(titles[j],fontsize=10)
fig.suptitle('TBM small geometry variants: fixed spherical envelope and all36 pipe ends\nBlue: coolant volume; orange: divider reference within coolant envelope. GEOMETRY ONLY, no CFD prediction',fontsize=12)
fig.tight_layout(rect=(0,0,1,.9));fig.savefig(O/'CAD-overview.png',dpi=180);plt.close(fig)

plane=cq.Face.makePlane(120,120,cq.Vector(0,0,20),cq.Vector(0,1,0))
base_water=water['B00_baseline'];base_metal=metal['B00_baseline']
def fill_section(ax,s,color):
 a=s.intersect(plane)
 if not a.Faces():return
 t=mesh(a,.006);ax.add_collection(PolyCollection(t[:,:,[0,2]],facecolor=color,edgecolor='none'))
def outline(ax,s,color,lw=.7):
 a=s.intersect(plane)
 for e in a.Edges():
  n=max(2,min(180,int(e.Length()/.012)+1));v=np.array([e.positionAt(t).toTuple() for t in np.linspace(0,1,n)])
  ax.plot(v[:,0],v[:,2],c=color,lw=lw)
fig,axs=plt.subplots(1,2,figsize=(9,6))
for ax,name,radius in zip(axs,['A01_root_R025','A02_root_R050'],[.25,.5]):
 fill_section(ax,base_water,'#b9deef');fill_section(ax,base_metal,'#9f9f9f')
 addition=base_water.cut(water[name]);fill_section(ax,addition,'#ec8b20');outline(ax,water[name],'#1d526e')
 ax.set(xlim=(21.2,25.2),ylim=(36.9,43),xlabel='X [mm]',ylabel='Z [mm]',title=f'Y=0 section: divider-root R{radius:.2f} mm');ax.set_aspect('equal');ax.grid(alpha=.15)
fig.legend(handles=[Patch(color='#b9deef',label='Water'),Patch(color='#9f9f9f',label='Original divider'),Patch(color='#ec8b20',label='Added divider material / removed water')],loc='lower center',ncol=3,fontsize=8)
fig.suptitle('Small fillets at the divider / spherical-wall junctions\nSame1 mm minimum divider thickness; no opening between water regions. GEOMETRY ONLY',fontsize=11)
fig.tight_layout(rect=(0,.07,1,.9));fig.savefig(O/'A-root-fillet-section.png',dpi=200);plt.close(fig)

cs=[sum(D['original_partition_face_x_mm'][k:k+2])/2 for k in [0,2]];c=cs[1];yy=np.linspace(-40,40,201)
fig,axs=plt.subplots(1,3,figsize=(13,4))
for ax,model,title in zip(axs,models[3:],['B01: inward0.50 mm','C01: right +1 deg','C02: right -1 deg']):
 par=model['parameters'];angle=par.get('right_partition_rotation_about_Z_deg',0);offset=par.get('inward_offset_each_mm',0)
 x=c-offset-np.tan(np.deg2rad(angle))*yy;half=.5/np.cos(np.deg2rad(angle))
 ax.fill_between(yy,c-.5,c+.5,color='#aaa',alpha=.45,label='Original1 mm divider')
 ax.fill_between(yy,x-half,x+half,color='#e98524',alpha=.5,label='Modified divider')
 ax.plot(yy,x,c='#a65b16',lw=1);ax.axvline(0,c='gray',ls=':',lw=.6)
 ax.set(xlim=(-40,40),ylim=(21.8,24.7),xlabel='Y [mm]: inlet side <0, outlet side >0',ylabel='X [mm]',title=title);ax.grid(alpha=.15)
axs[0].legend(fontsize=7)
fig.suptitle('XY projection of RIGHT divider; left divider is mirrored. Axes deliberately use different scales\nFixed pipes and spherical walls; tilted divider thickness remains1 mm measured normal to its plane',fontsize=11)
fig.tight_layout(rect=(0,0,1,.85));fig.savefig(O/'B-C-divider-position-angle.png',dpi=180);plt.close(fig)

import csv
ports=list(csv.DictReader((O/'ports.csv').open()))
fig,ax=plt.subplots(figsize=(7,7));theta=np.linspace(0,2*np.pi,400);ax.plot(47.9*np.cos(theta),47.9*np.sin(theta),c='gray',lw=.8)
for side,color in [('IN','#2874ad'),('OUT','#c54032')]:
 pp=[p for p in ports if p['side']==side];ax.scatter([float(p['x_mm']) for p in pp],[float(p['y_mm']) for p in pp],s=35,c=color,label=side)
 for p in pp:
  x=float(p['x_mm']);y=float(p['y_mm']);ax.text(x*1.13,y*1.13,p['port'],ha='center',va='center',fontsize=6,rotation=0)
for x0 in cs:ax.plot([x0,x0],[-41.5,41.5],c='#888',ls=':',lw=.7)
ax.set(xlim=(-60,60),ylim=(-60,60),xlabel='X [mm]',ylabel='Y [mm]',title='Unchanged36 port centres: CAD lookup labels\nVerify matching with existing STAR boundary names by position')
ax.set_aspect('equal');ax.legend(loc='upper right');ax.grid(alpha=.15);fig.tight_layout();fig.savefig(O/'unchanged-ports.png',dpi=180);plt.close(fig)
print('Rendered geometry-only overview, fillet section, divider-layout and port figures')
