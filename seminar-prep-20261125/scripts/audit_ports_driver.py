from pathlib import Path
import subprocess,sys,time
R=Path(__file__).resolve().parent
for name in ['B00','C01','C02']:
 while not (R/name/'mesh-end-UTC.txt').exists():time.sleep(5)
 p=subprocess.run([sys.executable,str(R/'audit_ports.py'),name],cwd=R,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (R/f'log.port-audit-{name}').write_text(p.stdout)
 if p.returncode:raise RuntimeError(name+' port audit failed; inspect log')
 print(name+' native port geometry audited',flush=True)
(R/'PORT_AUDITS_COMPLETE.txt').write_text('All three pure-geometry port audits completed; no field arrays or flow solver.\n')
