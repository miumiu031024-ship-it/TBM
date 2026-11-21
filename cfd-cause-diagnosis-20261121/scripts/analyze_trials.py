"""One-factor steady tests. Iterations are not physical time."""
from pathlib import Path
import re,json,csv,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;O=R/'diagnosis';O.mkdir(exist_ok=True)
rho=746.200955901462
names=json.loads((R/'trials.json').read_text())
def read_u(path):
 b=path.read_bytes();m=re.search(rb'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(',b);assert m
 return np.frombuffer(b,dtype='<f8',offset=m.end(),count=3*int(m.group(1))).reshape(-1,3)
def read_probes(path,vec=False):
 if not path.exists():return {}
 rows={}
 for line in path.read_text().splitlines():
  if line.startswith('#') or not line.strip():continue
  a=np.fromstring(line.replace('(',' ').replace(')',' '),sep=' ')
  if len(a)==1+1359*(3 if vec else 1):rows[a[0]]=a[1:].reshape(-1,3) if vec else a[1:]
 return rows
baseu=read_u(R/'SeedSide/500/U');summaries={};allhistory=[]
fig,axs=plt.subplots(2,2,figsize=(12,7))
for name in names:
 C=R/name;lf=C/'log.simpleFoam'
 if not lf.exists():continue
 log=lf.read_text();res=[]
 for b in re.split(r'\nTime = ',log)[1:]:
  t=float(b.splitlines()[0]);vv={}
  for v,x in re.findall(r'Solving for (\w+), Initial residual = ([\deE+.-]+)',b):vv.setdefault(v,float(x))
  if all(k in vv for k in ['Ux','Uy','Uz','p','k','omega']):res.append(dict(iteration=t,**vv))
 if not res:continue
 with (O/f'{name}-residuals.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(res[0]));w.writeheader();w.writerows(res)
 mon={}
 for side in ['inlet','outlet']:
  for i in range(1,7):
   for suffix in ['Flux','Pressure']:
    key=f'{side}{i:02d}{suffix}';rows={}
    for p in (C/'postProcessing'/key).glob('*/surfaceFieldValue.dat'):
     for ln in p.read_text().splitlines():
      if ln and not ln.startswith('#'):
       a=list(map(float,ln.split()));rows[a[0]]=a[1]
    mon[key]=rows
 its=sorted(set.intersection(*[set(v) for v in mon.values()])) if all(mon.values()) else []
 hist=[]
 for t in its:
  qo=np.array([mon[f'outlet{i:02d}Flux'][t]*rho for i in range(1,7)]);qi=-sum(mon[f'inlet{i:02d}Flux'][t]*rho for i in range(1,7))
  hist.append(dict(case=name,iteration=t,inlet_kg_s=qi,outlet_kg_s=float(qo.sum()),outlet_spread_percent=float(np.ptp(qo)/qo.mean()*100),inlet_unweighted_static_Pa=float(np.mean([mon[f'inlet{i:02d}Pressure'][t]*rho for i in range(1,7)])),**{f'outlet{i+1:02d}_kg_s':float(q) for i,q in enumerate(qo)}))
 allhistory.extend(hist)
 pu=read_probes(C/'postProcessing/interiorProbes/500/U',True)
 pv=read_probes(C/'postProcessing/interiorProbes/500/vorticity',True)
 pp=read_probes(C/'postProcessing/interiorProbes/500/p')
 pr=[];prev=None;first=None
 for t,u in sorted(pu.items()):
  valid=np.isfinite(u).all(axis=1)&(abs(u)<1e10).all(axis=1)
  if prev is None:prev=(t,u,valid);first=(t,u,valid)
  good=valid&prev[2];speed=np.linalg.norm(u[valid],axis=1)
  change=np.sqrt(((u[good]-prev[1][good])**2).sum()/(prev[1][good]**2).sum())
  vort=pv.get(t);pressure=pp.get(t)
  row=dict(case=name,iteration=t,valid_points=int(valid.sum()),mean_speed_m_s=float(speed.mean()),low_fraction_0p5=float((speed<.5).mean()),vector_L2_since_previous=float(change))
  if vort is not None:
   mask=valid&(abs(vort)<1e10).all(axis=1);row['mean_vorticity_s_inv']=float(np.linalg.norm(vort[mask],axis=1).mean())
  if pressure is not None:row['probe_mean_pressure_Pa']=float(pressure[valid].mean()*rho)
  pr.append(row);prev=(t,u,valid)
 if pr:
  keys=list(dict.fromkeys(k for r in pr for k in r))
  with (O/f'{name}-internal-monitor.csv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(pr)
 result=dict(last_iteration=res[-1]['iteration'],completed_normally='\nEnd\n' in log,final_residuals=res[-1],residual_goals_met=all(res[-1][k]<(1e-4 if k=='p' else 1e-5) for k in ['p','Ux','Uy','Uz','k','omega']),final_boundary_monitors=hist[-1] if hist else None)
 if (C/'600/U').exists():
  u=read_u(C/'600/U');up=read_u(C/'550/U')
  result['native_velocity_relative_L2_500_to_600']=float(np.sqrt(((u-baseu)**2).sum()/(baseu**2).sum()))
  result['native_velocity_relative_L2_550_to_600']=float(np.sqrt(((u-up)**2).sum()/(up**2).sum()))
 if len(pr)>1:
  recent=[r for r in pr if r['iteration']>=550]
  result['recent_low_fraction_range_percentage_points']=100*float(np.ptp([r['low_fraction_0p5'] for r in recent]))
  result['maximum_recent_5iteration_vector_change']=max(r['vector_L2_since_previous'] for r in recent)
 summaries[name]=result
 t=[r['iteration']-500 for r in res];axs[0,0].semilogy(t,[r['p'] for r in res],label=name);axs[0,1].semilogy(t,[r['Ux'] for r in res],label=name)
 if hist:
  axs[1,0].plot([r['iteration']-500 for r in hist],[r['outlet_spread_percent'] for r in hist],label=name)
 if pr:axs[1,1].plot([r['iteration']-500 for r in pr],[100*r['low_fraction_0p5'] for r in pr],label=name)
for ax,label in zip(axs.ravel(),['Initial p residual','Initial Ux residual','6-outlet spread / mean [%]','Internal probes speed <0.5 m/s [%]']):
 ax.set(xlabel='Extra steady iterations from identical restart',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.suptitle('One-factor numerical tests: identical mesh, geometry, BCs, and restart field\nIteration histories are not physical time histories',fontsize=11);fig.tight_layout();fig.savefig(O/'one-factor-comparison.png',dpi=180);plt.close(fig)
(O/'steady-trials-summary.json').write_text(json.dumps(summaries,indent=2))
if allhistory:
 with (O/'steady-boundary-history.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(allhistory[0]));w.writeheader();w.writerows(allhistory)
print(json.dumps(summaries,indent=2))
