"""Localize measured steady-iteration updates; not a physical-time/cause inference."""
from pathlib import Path
import re,json,csv
import numpy as np
R=Path(__file__).resolve().parent;O=R/'results';rho=746.200955901462
rows=[]
for name in ['F01_upwind_long','F02_second_order_restore']:
 C=R/name
 if not C.exists():continue
 data={f:{} for f in ['U','p','vorticity']};xyz=None
 for field in data:
  for p in sorted((C/'postProcessing/interiorProbes').glob('*/'+field),key=lambda p:float(p.parent.name)):
   body=p.read_text()
   if xyz is None:xyz=np.array([list(map(float,s.split())) for s in re.findall(r'^# Probe \d+ \(([^)]+)\)',body,re.M)]);assert xyz.shape==(1359,3)
   for line in body.splitlines():
    if not line.strip() or line.startswith('#'):continue
    a=np.fromstring(line.replace('(',' ').replace(')',' '),sep=' ');vec=field!='p'
    if len(a)==1+1359*(3 if vec else 1):data[field][a[0]]=a[1:].reshape(-1,3) if vec else a[1:]
 if xyz is None:continue
 regions={'All':np.ones(1359,bool),'Lower_z_under_5mm':xyz[:,2]<.005,'Middle_z_5_to_20mm':(xyz[:,2]>=.005)&(xyz[:,2]<.020),'Upper_z_over_20mm':xyz[:,2]>=.020,'Within_1mm_of_divider':abs(xyz[:,0]+.0236999776773)<.001,'Inlet_half_y_negative':xyz[:,1]<0,'Outlet_half_y_positive':xyz[:,1]>0}
 for p in sorted((O/name).glob('assessment-*.json')):
  q=json.loads(p.read_text());end=q['iteration'];start=end-100
  if not all(t in data['U'] for t in [start,end]):continue
  for region,mask in regions.items():
   a,b=data['U'][start][mask],data['U'][end][mask];assert np.isfinite(a).all() and np.isfinite(b).all() and (abs(b)<1e10).all()
   diff=np.linalg.norm(b-a,axis=1);al=np.linalg.norm(a,axis=1)<.5;bl=np.linalg.norm(b,axis=1)<.5;union=al|bl
   row=dict(case=name,iteration=end,from_iteration=start,region=region,sampled_points=int(mask.sum()),mean_vector_difference_m_s=float(diff.mean()),relative_vector_L2_change=float(np.sqrt(((b-a)**2).sum()/(a**2).sum())),low_mask_IoU=float((al&bl).sum()/union.sum()) if union.any() else 1.)
   if end in data['vorticity']:row['mean_vorticity_magnitude_s_inv']=float(np.linalg.norm(data['vorticity'][end][mask],axis=1).mean())
   if start in data['p'] and end in data['p']:row['mean_sample_static_pressure_change_Pa']=float(rho*(data['p'][end][mask]-data['p'][start][mask]).mean())
   rows.append(row)
if rows:
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (O/'spatial-update-location.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 print('Saved regional steady-iteration update diagnostics',len(rows))
