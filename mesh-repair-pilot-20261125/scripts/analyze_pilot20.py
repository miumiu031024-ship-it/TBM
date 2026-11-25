from pathlib import Path
import re,json,csv,sys,gzip,hashlib,shutil
import numpy as np
R=Path(__file__).resolve().parent;O=Path('/workspace/TBM/mesh-repair-pilot-20261125');rho=746.200955901462
n=sys.argv[1];C=R/(n+'_R07');D=O/n;D.mkdir(exist_ok=True,parents=True)
def savecsv(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
log=(C/'log.simpleFoam20').read_text();assert re.findall(r'^Time = (\d+)',log,re.M)==[str(i) for i in range(1,21)];assert '\nEnd\n' in log
res=[]
for b in re.split(r'\nTime = ',log)[1:]:
 i=int(b.splitlines()[0]);row=dict(iteration=i)
 for v,x,y,z in re.findall(r'Solving for (\w+), Initial residual = ([\deE+.-]+), Final residual = ([\deE+.-]+), No Iterations (\d+)',b):
  row.setdefault(v,float(x));row.setdefault(v+'_linear_final',float(y));row.setdefault(v+'_linear_iterations',int(z))
 assert all(v in row for v in ['p','Ux','Uy','Uz','k','omega']);res.append(row)
savecsv(D/'residuals20.csv',res)
mon={};areas={}
for side in ['inlet','outlet']:
 for j in range(1,19):
  for suffix in ['Flux','Pressure']:
   key=f'{side}{j:02d}{suffix}';rows={}
   for p in sorted((C/'postProcessing'/key).glob('*/surfaceFieldValue.dat')):
    body=p.read_text();m=re.search(r'#\s*Area\s*:\s*([\deE+.-]+)',body)
    if m:areas[f'{side}{j:02d}']=float(m.group(1))
    for line in body.splitlines():
     if line and not line.startswith('#'):
      a=list(map(float,line.split()));rows[int(a[0])]=a[1]
   assert 20 in rows,key;mon[key]=rows
hist=[];outrows=[];regions=[]
for i in range(1,21):
 qi=np.array([-mon[f'inlet{j:02d}Flux'][i]*rho for j in range(1,19)]);qo=np.array([mon[f'outlet{j:02d}Flux'][i]*rho for j in range(1,19)])
 pin=rho*np.average([mon[f'inlet{j:02d}Pressure'][i] for j in range(1,19)],weights=[areas[f'inlet{j:02d}'] for j in range(1,19)])
 pout=rho*np.average([mon[f'outlet{j:02d}Pressure'][i] for j in range(1,19)],weights=[areas[f'outlet{j:02d}'] for j in range(1,19)])
 row=dict(iteration=i,inlet_kg_s=float(qi.sum()),outlet_kg_s=float(qo.sum()),net_imbalance_kg_s=float(qo.sum()-qi.sum()),net_imbalance_percent=float(100*(qo.sum()-qi.sum())/qi.sum()),inlet_mean_static_Pa=float(pin),outlet_mean_static_Pa=float(pout),static_pressure_difference_Pa=float(pin-pout))
 for k,label in enumerate(['Water_Xminus','Water_Center','Water_Xplus']):
  s=slice(6*k,6*k+6);row[label+'_imbalance_percent']=float(100*(qo[s].sum()-qi[s].sum())/qi[s].sum());row[label+'_spread_percent']=float(100*np.ptp(qo[s])/qo[s].mean())
  if i==20:regions.append(dict(region=label,inlet_kg_s=float(qi[s].sum()),outlet_kg_s=float(qo[s].sum()),imbalance_percent=row[label+'_imbalance_percent'],outlet_spread_percent=row[label+'_spread_percent'],first_outlet_kg_s=float(qo[6*k]),sixth_outlet_kg_s=float(qo[6*k+5])))
 row.update({f'outlet{j:02d}_kg_s':float(q) for j,q in enumerate(qo,1)});hist.append(row)
 if i==20:outrows=[dict(outlet=f'OUT{j:02d}',mass_flow_kg_s=float(q),NOTE='20 iteration startup only; not converged outlet distribution') for j,q in enumerate(qo,1)]
savecsv(D/'boundary-history20.csv',hist);savecsv(D/'18-outlets20.csv',outrows);savecsv(D/'regional-flow20.csv',regions)
fields={}
for f in ['U','p','k','omega','nut','phi']:
 p=C/'20'/(f+'.gz')
 if not p.exists():
  raw=C/'20'/f;assert raw.exists(),raw
  with raw.open('rb') as src,gzip.open(p,'wb',compresslevel=1) as dst:shutil.copyfileobj(src,dst)
  # Check gzip round trip hash before removing only new duplicate uncompressed checkpoint.
  h=hashlib.sha256()
  with gzip.open(p,'rb') as z:
   for block in iter(lambda:z.read(8*1024*1024),b''):h.update(block)
  assert h.hexdigest()==hashlib.sha256(raw.read_bytes()).hexdigest();raw.unlink()
 assert p.exists(),p
 with gzip.open(p,'rb') as g:b=g.read()
 m=re.search(rb'internalField\s+nonuniform\s+List<(scalar|vector)>\s+(\d+)\s*\(',b);assert m,(n,f)
 count=int(m.group(2));comp=3 if m.group(1)==b'vector' else 1;a=np.frombuffer(b,dtype='<f8',offset=m.end(),count=count*comp)
 assert np.isfinite(a).all(),(n,f,'nonfinite');q=dict(internal_entities=count,min=float(a.min()),max=float(a.max()),all_finite=True,compressed_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 if f=='U':q.update(max_speed_m_s=float(np.linalg.norm(a.reshape(-1,3),axis=1).max()))
 if f in ['k','omega']:assert q['min']>0,(n,f,q)
 if f=='nut':assert q['min']>=0
 fields[f]=q
peakrows=list(csv.DictReader((R/'pilot-memory-samples.csv').open()));samples=[q for q in peakrows if q['case']==n+'_R07' and float(q['solver_sum_RSS_KiB'])>0]
assert samples,'No solver memory samples'
assert all(int(q['mesh_check_process_count'])==0 for q in samples),'Concurrent mesh checker polluted benchmark'
result=dict(model=n,iterations=20,normal_end=True,all_fields_finite=True,fields=fields,final_initial_residuals={v:res[-1][v] for v in ['p','Ux','Uy','Uz','k','omega']},start_boundary=hist[0],last_boundary=hist[-1],regions=regions,
 startup_yPlus_walls=[float(v) for v in re.search(r'patch walls y\+ : min = ([\deE+.-]+), max = ([\deE+.-]+), average = ([\deE+.-]+)',log).groups()],
 solver_elapsed_seconds=float(re.findall(r'ClockTime = ([\deE+.-]+) s',log)[-1]),driver_time=(C/'solver-driver-time.txt').read_text().strip(),
 sampled_peak_solver_summed_RSS_GiB=max(float(q['solver_sum_RSS_KiB'])/1024**2 for q in samples),sampled_peak_container_memory_including_cache_GiB=max(float(q['container_memory_current_bytes'])/1024**3 for q in samples),sampled_peak_container_anon_GiB=max(float(q['container_anon_bytes'])/1024**3 for q in samples),memory_sampling_interval_seconds=5,
 warnings=['\n'.join(log.splitlines()[i:i+6]) for i,line in enumerate(log.splitlines()) if '--> FOAM Warning' in line or 'bounding ' in line or 'FOAM FATAL' in line],actual_solver_fatal=bool(re.search(r'FOAM FATAL|\bSIGFPE\b.*[Ss]ignal|[Ff]loating point exception \(core',log)),
 last10_max_abs_global_imbalance_percent=max(abs(q['net_imbalance_percent']) for q in hist[-10:]),last10_inlet_static_pressure_range_Pa=max(q['inlet_mean_static_Pa'] for q in hist[-10:])-min(q['inlet_mean_static_Pa'] for q in hist[-10:]),
 NOTE='Startup test only. Steady iteration is NOT physical time. No design improvement or convergence established. 200 iterations NOT executed.')
meta=json.loads((C/'case-metadata.json').read_text());meta.update(flow_solver_executed=True,completed_steady_iterations=20,solver='simpleFoam',normal_solver_end=True);(C/'case-metadata.json').write_text(json.dumps(meta,indent=2))
(D/'startup-summary.json').write_text(json.dumps(result,indent=2));print(n,'20 iterations finite, all six fields preserved, runtime',result['solver_elapsed_seconds'],'s, peak summed RSS',result['sampled_peak_solver_summed_RSS_GiB'],'GiB',flush=True)
# Reconstructed global state is preserved. Only this task\'s newly generated redundant decomposition is removed.
procfiles={str(p.relative_to(C)):hashlib.sha256(p.read_bytes()).hexdigest() for d in C.glob('processor[0-9]*') for p in d.rglob('*') if p.is_file() and p.parent.name=='20'}
(C/'reconstructed-processor20-sha256.json').write_text(json.dumps(procfiles,indent=2))
assert len(procfiles)>=24
for d in C.glob('processor[0-9]*'):
 assert d.is_dir() and not d.is_symlink();shutil.rmtree(d)
(C/'REDUNDANT_DECOMPOSITION_REMOVED.txt').write_text('Only new processor copies removed AFTER all reconstructed global time20 fields verified finite, positive turbulence quantities, and SHA256 recorded. Native global mesh, 0, time20 fields, postProcessing and logs retained. Originals never modified.\n')
