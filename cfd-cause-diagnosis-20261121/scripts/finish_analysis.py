"""Postprocess both normally completed pilots; do not package or publish."""
from pathlib import Path
import subprocess,time,os,sys
R=Path(__file__).resolve().parent
for i in range(180):
 f=R/'log.transient-driver';s=f.read_text() if f.exists() else ''
 if 'T40_dt1e5 complete' in s:break
 if f.exists() and ('FOAM FATAL' in s or 'exited on signal' in s):raise RuntimeError('Transient driver failed')
 time.sleep(10)
else:raise TimeoutError('Pilot completion timeout')
for n in ['T30_dt2e5','T40_dt1e5']:assert '\nEnd\n' in (R/n/'log.pimpleFoam').read_text()
env=os.environ|{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VTK_SMP_MAX_THREADS':'1'}
for script,args in [('analyze_transient.py',[]),('diagnose_sections.py',['T40_dt1e5']),('diagnose_flux_consistency.py',[]),('export_physical_samples.py',[]),('write_report.py',[])]:
 with (R/('log.final-'+script+'.txt')).open('w') as f:subprocess.run([sys.executable,script,*args],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 print('Completed',script,flush=True)
print('Ready for figure/report review; not packaged or published',flush=True)
