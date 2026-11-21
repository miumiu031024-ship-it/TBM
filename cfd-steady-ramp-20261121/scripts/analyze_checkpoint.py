"""Assess nonlinear steady iteration convergence using independent observables."""
from pathlib import Path
import json,re,csv,sys,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-ramp-plots')
import numpy as np
R=Path(__file__).resolve().parent

def read_u(path):
 b=path.read_bytes();m=re.search(rb'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(',b);assert m
 a=np.frombuffer(b,dtype='<f8',offset=m.end(),count=3*int(m.group(1))).reshape(-1,3);assert np.isfinite(a).all();return a

def csv_write(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def assess(name,end):
 C=R/name;O=R/'results'/name;O.mkdir(parents=True,exist_ok=True)
 plan=json.loads((R/'plan.json').read_text())['assessment'];rho=746.200955901462;res={}
 for p in C.glob('log.simpleFoam-*'):
  for block in re.split(r'\nTime = ',p.read_text())[1:]:
   t=float(block.splitlines()[0]);vals={}
   for v,x in re.findall(r'Solving for (\w+), Initial residual = ([\deE+.-]+)',block):vals.setdefault(v,float(x))
   if all(k in vals for k in ['p','Ux','Uy','Uz','k','omega']):res[t]=dict(iteration=t,**{k:vals[k] for k in ['p','Ux','Uy','Uz','k','omega']})
 assert end in res
 res=[res[t] for t in sorted(res) if t<=end];csv_write(O/'residuals.csv',res)
 recent=[q for q in res if q['iteration']>end-100];assert len(recent)==100
 residual_max={k:max(q[k] for q in recent) for k in plan['residual_max_last100']}
 mon={};areas={}
 for side in ['inlet','outlet']:
  for i in range(1,7):
   for suffix in ['Flux','Pressure']:
    key=f'{side}{i:02d}{suffix}';rows={}
    for p in sorted((C/'postProcessing'/key).glob('*/surfaceFieldValue.dat'),key=lambda p:float(p.parent.name)):
     body=p.read_text();areas[f'{side}{i:02d}']=float(re.search(r'#\s*Area\s*:\s*([\deE+.-]+)',body).group(1))
     for ln in body.splitlines():
      if ln and not ln.startswith('#'):
       a=list(map(float,ln.split()));rows[a[0]]=a[1]
    assert rows,key;mon[key]=rows
 times=sorted(set.intersection(*[set(q) for q in mon.values()]));times=[t for t in times if t<=end];hist=[]
 for t in times:
  qo=np.array([mon[f'outlet{i:02d}Flux'][t]*rho for i in range(1,7)]);qi=-sum(mon[f'inlet{i:02d}Flux'][t]*rho for i in range(1,7))
  pin=rho*np.average([mon[f'inlet{i:02d}Pressure'][t] for i in range(1,7)],weights=[areas[f'inlet{i:02d}'] for i in range(1,7)])
  hist.append(dict(iteration=t,inlet_kg_s=qi,outlet_kg_s=float(qo.sum()),mass_imbalance_fraction=float((qi-qo.sum())/qi),outlet_spread_percent=float(np.ptp(qo)/qo.mean()*100),inlet_area_mean_static_Pa=float(pin),**{f'outlet{i+1:02d}_kg_s':float(q) for i,q in enumerate(qo)}))
 csv_write(O/'boundary-history.csv',hist);h=[q for q in hist if q['iteration']>=end-100];assert len(h)>=20 and h[-1]['iteration']==end
 flow_ranges={f'outlet{i:02d}':float(np.ptp([q[f'outlet{i:02d}_kg_s'] for q in h])/np.mean([q[f'outlet{i:02d}_kg_s'] for q in h])) for i in range(1,7)}
 probes={}
 for p in sorted((C/'postProcessing/interiorProbes').glob('*/U'),key=lambda p:float(p.parent.name)):
  for line in p.read_text().splitlines():
   if not line.strip() or line.startswith('#'):continue
   a=np.fromstring(line.replace('(',' ').replace(')',' '),sep=' ')
   if len(a)==1+1359*3 and a[0]<=end:probes[a[0]]=a[1:].reshape(-1,3)
 pt=sorted(probes);pr=[];prev=None
 for t in pt:
  u=probes[t];valid=np.isfinite(u).all(axis=1)&(abs(u)<1e10).all(axis=1);assert valid.sum()==1359
  speed=np.linalg.norm(u,axis=1);lo=speed<.5
  row=dict(iteration=t,valid_points=1359,mean_speed_m_s=float(speed.mean()),low_fraction_0p5=float(lo.mean()),vector_L2_since_previous=None,low_mask_IoU_since_previous=None)
  if prev is not None:
   old=prev[1];ol=np.linalg.norm(old,axis=1)<.5;union=lo|ol
   row['vector_L2_since_previous']=float(np.sqrt(((u-old)**2).sum()/(old**2).sum()))
   row['low_mask_IoU_since_previous']=float((lo&ol).sum()/union.sum()) if union.any() else 1.
  pr.append(row);prev=(t,u)
 csv_write(O/'internal-monitor.csv',pr)
 recent_pt=[t for t in pt if t>=end-100];assert len(recent_pt)>=6
 changes=[]
 for ta,tb in zip(recent_pt[:-1],recent_pt[1:]):
  a,b=probes[ta],probes[tb];changes.append(float(np.sqrt(((b-a)**2).sum()/(a**2).sum())))
 a,b=probes[recent_pt[0]],probes[recent_pt[-1]];al=np.linalg.norm(a,axis=1)<.5;bl=np.linalg.norm(b,axis=1)<.5;un=al|bl;iou=float((al&bl).sum()/un.sum()) if un.any() else 1.
 low_range=100*float(np.ptp([q['low_fraction_0p5'] for q in pr if q['iteration']>=end-100]))
 ua=read_u(C/str(end-100)/'U');ub=read_u(C/str(end)/'U');assert ua.shape==ub.shape==(1264925,3)
 native_change=float(np.sqrt(((ub-ua)**2).sum()/(ua**2).sum()))
 gates=dict(residuals=all(residual_max[k]<goal for k,goal in plan['residual_max_last100'].items()),native_velocity=native_change<plan['native_velocity_relative_L2_last100_below'],probe_velocity=max(changes)<plan['maximum_20_iteration_probe_vector_relative_L2_last100_below'],low_fraction=low_range<plan['probe_low_fraction_range_last100_pp_below'],low_location=iou>=plan['low_mask_IoU_last100_at_least'],outlet_flows=max(flow_ranges.values())<plan['each_outlet_flow_range_last100_fraction_below'],mass_balance=max(abs(q['mass_imbalance_fraction']) for q in h)<plan['absolute_total_mass_imbalance_fraction_below'])
 summary=dict(case=name,iteration=end,normal_solver_end=all('\nEnd\n' in p.read_text() for p in C.glob('log.simpleFoam-*')),last_residuals=res[-1],maximum_residuals_last100=residual_max,native_velocity_L2_last100=native_change,maximum_20iter_probe_velocity_L2_last100=max(changes),low_fraction_range_last100_pp=low_range,low_mask_IoU_last100=iou,each_outlet_flow_range_fraction_last100=flow_ranges,final_monitors=hist[-1],criteria=gates,all_steady_numerical_criteria_met=all(gates.values()),grid_independence_established=False,physical_time_interpretation=False)
 (O/f'assessment-{end}.json').write_text(json.dumps(summary,indent=2));(O/'latest-assessment.json').write_text(json.dumps(summary,indent=2))
 print(json.dumps(summary,indent=2),flush=True);return summary
if __name__=='__main__':assess(sys.argv[1],int(sys.argv[2]))
