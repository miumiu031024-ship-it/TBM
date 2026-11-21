"""Wait for normally completed/exported stages, then analyze. No publication."""
from pathlib import Path
import subprocess,time,os,sys
R=Path(__file__).resolve().parent
for i in range(1080):
 if (R/'SOLVES_COMPLETE.json').exists():break
 p=R/'log.driver'
 if p.exists() and ('Traceback (most recent call last)' in p.read_text() or 'Floating point exception (8)' in p.read_text()):raise RuntimeError('Driver failed; inspect log.driver')
 time.sleep(10)
else:raise TimeoutError('Stage completion timeout')
env=os.environ|{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VTK_SMP_MAX_THREADS':'1'}
for script in ['final_analysis.py','export_monitor_history.py','spatial_monitor_analysis.py','write_report.py']:
 with (R/('log.final-'+script)).open('w') as f:subprocess.run([sys.executable,script],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 print('Completed',script,flush=True)
print('Ready for human-readable output review; not packaged or published',flush=True)
