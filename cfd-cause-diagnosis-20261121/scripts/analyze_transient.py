"""Physical-time pilot diagnostics, explicitly separate from steady iterations."""
from pathlib import Path
import json,re,csv,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-cause-plots')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;O=R/'diagnosis';rho=746.200955901462
names=['T30_dt2e5','T40_dt1e5']
def read_u(path):
 b=path.read_bytes();m=re.search(rb'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(',b);assert m
 return np.frombuffer(b,dtype='<f8',offset=m.end(),count=3*int(m.group(1))).reshape(-1,3)
def probe_values(path,vec=False):
 if not path.exists():return {}
 rows={}
 for line in path.read_text().splitlines():
  if line.startswith('#') or not line.strip():continue
  a=np.fromstring(line.replace('(',' ').replace(')',' '),sep=' ')
  if len(a)==1+1359*(3 if vec else 1):rows[round(a[0],10)]=a[1:].reshape(-1,3) if vec else a[1:]
 return rows
summaries={};allprobe=[];allflows=[]
fig,axs=plt.subplots(2,2,figsize=(12,7))
for name in names:
 C=R/name;f=C/'log.pimpleFoam'
 if not f.exists():continue
 log=f.read_text();records=[];co=None;row=None
 for line in log.splitlines():
  cm=re.match(r'Courant Number mean: ([\deE+.-]+) max: ([\deE+.-]+)',line)
  if cm:co=tuple(map(float,cm.groups()))
  tm=re.match(r'Time = (\S+)',line)
  if tm:
   if row and 'p_last_final' in row:records.append(row)
   row=dict(time_s=float(tm.group(1)),courant_mean=co[0] if co else None,courant_max=co[1] if co else None)
  sm=re.search(r'Solving for (\w+), Initial residual = ([\deE+.-]+), Final residual = ([\deE+.-]+)',line)
  if sm and row is not None:
   v,ini,fin=sm.groups();row.setdefault(v+'_first_initial',float(ini));row[v+'_last_initial']=float(ini);row[v+'_last_final']=float(fin)
  cm=re.search(r'continuity errors : sum local = ([\deE+.-]+), global = ([\deE+.-]+)',line)
  if cm and row is not None:row['continuity_local_last']=float(cm.group(1));row['continuity_global_last']=float(cm.group(2))
 if row and 'p_last_final' in row:records.append(row)
 if not records:continue
 keys=list(dict.fromkeys(k for q in records for k in q))
 with (O/f'{name}-step-diagnostics.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(records)
 mon={}
 for side in ['inlet','outlet']:
  for i in range(1,7):
   key=f'{side}{i:02d}';vals={}
   for p in (C/'postProcessing'/(key+'Flux')).glob('*/surfaceFieldValue.dat'):
    for line in p.read_text().splitlines():
     if line and not line.startswith('#'):
      a=list(map(float,line.split()));vals[round(a[0],10)]=a[1]*rho
   mon[key]=vals
 times=sorted(set.intersection(*[set(q) for q in mon.values()])) if all(mon.values()) else []
 flows=[]
 for t in times:
  qi=-sum(mon[f'inlet{i:02d}'][t] for i in range(1,7));qs=np.array([mon[f'outlet{i:02d}'][t] for i in range(1,7)])
  flows.append(dict(case=name,time_s=t,inlet_kg_s=qi,outlet_kg_s=float(qs.sum()),mass_imbalance_fraction=float((qi-qs.sum())/qi),flow_spread_percent=float(np.ptp(qs)/qs.mean()*100),outlet06_minus_01_kg_s=float(qs[-1]-qs[0]),**{f'outlet{i+1:02d}_kg_s':float(q) for i,q in enumerate(qs)}))
 allflows.extend(flows)
 pu=probe_values(C/'postProcessing/interiorProbes/0/U',True);pv=probe_values(C/'postProcessing/interiorProbes/0/vorticity',True);pp=probe_values(C/'postProcessing/interiorProbes/0/p')
 probes=[];prev=None
 for t,u in sorted(pu.items()):
  good=np.isfinite(u).all(axis=1)&(abs(u)<1e10).all(axis=1)
  row=dict(case=name,time_s=t,valid_points=int(good.sum()),low_fraction_0p5=float((np.linalg.norm(u[good],axis=1)<.5).mean()),mean_speed_m_s=float(np.linalg.norm(u[good],axis=1).mean()))
  if prev is not None:
   g=good&prev[1];row['vector_relative_L2_since_previous']=float(np.sqrt(((u[g]-prev[0][g])**2).sum()/(prev[0][g]**2).sum()))
  if t in pv:row['mean_vorticity_s_inv']=float(np.linalg.norm(pv[t][good],axis=1).mean())
  if t in pp:row['mean_sample_static_Pa']=float(pp[t][good].mean()*rho)
  probes.append(row);prev=(u,good)
 allprobe.extend(probes)
 dt=json.loads((C/'case-metadata.json').read_text())['delta_t_s']
 summaries[name]=dict(completed_normally='\nEnd\n' in log,last_physical_time_s=records[-1]['time_s'],delta_t_s=dt,completed_steps=len(records),maximum_courant_number=max(q['courant_max'] for q in records),final_step=records[-1],final_flows=flows[-1] if flows else None,physical_probe_frames=len(probes),mean_residence_time_estimate_s=rho*1.0776389e-5/(.23/3),periodicity_established=False,statistically_converged_mean_established=False)
 t=[q['time_s'] for q in records];axs[0,0].plot(t,[q['courant_max'] for q in records],label=name)
 if flows:
  axs[0,1].plot([q['time_s'] for q in flows],[q['flow_spread_percent'] for q in flows],label=name)
 if probes:
  axs[1,0].plot([q['time_s'] for q in probes],[q['low_fraction_0p5']*100 for q in probes],label=name)
  axs[1,1].plot([q['time_s'] for q in probes],[q.get('mean_vorticity_s_inv',np.nan) for q in probes],label=name)
for ax,y in zip(axs.ravel(),['Maximum Courant number','6-outlet spread / mean [%]','Probe fraction speed <0.5 m/s [%]','Probe mean vorticity magnitude [1/s]']):ax.set(xlabel='Physical time [s]',ylabel=y);ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.suptitle('Two time-step URANS pilots from the same unsettled initial field\nStartup interval only: not statistically converged means or proof of a physical period',fontsize=11);fig.tight_layout();fig.savefig(O/'physical-time-pilot-comparison.png',dpi=180);plt.close(fig)
for file,rows in [('physical-probe-history.csv',allprobe),('physical-boundary-history.csv',allflows)]:
 if rows:
  keys=list(dict.fromkeys(k for r in rows for k in r))
  with (O/file).open('w') as f:
   w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
if all((R/n/'.002/U').exists() or (R/n/'0.002/U').exists() for n in names):
 a,b=[read_u(R/n/'0.002/U') for n in names]
 summaries['paired_final_velocity_relative_L2']=float(np.sqrt(((a-b)**2).sum()/(b**2).sum()))
if all(n in summaries and summaries[n]['completed_normally'] for n in names):
 ua,ub=[probe_values(R/n/'postProcessing/interiorProbes/0/U',True)[.002] for n in names]
 good=np.isfinite(ua).all(axis=1)&np.isfinite(ub).all(axis=1)&(abs(ua)<1e10).all(axis=1)&(abs(ub)<1e10).all(axis=1)
 summaries['paired_final_probe_differences']=dict(valid_points=int(good.sum()),vector_relative_L2=float(np.sqrt(((ua[good]-ub[good])**2).sum()/(ub[good]**2).sum())),maximum_vector_difference_m_s=float(np.linalg.norm(ua[good]-ub[good],axis=1).max()),low_fraction_0p5_change_percentage_points=float(100*((np.linalg.norm(ua[good],axis=1)<.5).mean()-(np.linalg.norm(ub[good],axis=1)<.5).mean())))
 summaries['measured_clock_seconds']={n:float(re.findall(r'ClockTime = ([\deE+.-]+) s',(R/n/'log.pimpleFoam').read_text())[-1]) for n in names}
(O/'transient-pilots-summary.json').write_text(json.dumps(summaries,indent=2))
print(json.dumps(summaries,indent=2))
