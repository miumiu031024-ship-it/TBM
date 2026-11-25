from pathlib import Path
import time,subprocess,sys,datetime
R=Path(__file__).resolve().parent
while not ((R/'MESHES_COMPLETE.txt').exists() and (R/'PORT_AUDITS_COMPLETE.txt').exists() and 'Memory sampling complete' in (R/'log.memory-monitor').read_text()):time.sleep(5)
for script in ['inspect_quality_cells.py','audit_regions.py','audit_meshes.py','plot_native_ports.py','write_final_report.py','make_slides.py']:
 with (R/f'log.final-{script}').open('w') as f:
  p=subprocess.run([sys.executable,str(R/script)],stdout=f,stderr=subprocess.STDOUT,cwd=R)
 if p.returncode:raise RuntimeError('Failed '+script+'; inspect corresponding log')
 print('Completed '+script,flush=True)
(R/'PREPARATION_ANALYSIS_COMPLETE.txt').write_text(datetime.datetime.now(datetime.timezone.utc).isoformat()+'\nNo flow solve; output ready for visual review, NOT yet packaged or published.\n')
print('Prepared outputs ready for review; not packaged/published',flush=True)
