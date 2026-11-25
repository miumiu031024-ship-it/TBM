from pathlib import Path
import time,subprocess,sys,json,datetime
R=Path(__file__).resolve().parent
for n in ['B00','C01','C02']:
 while not (R/(n+'_R07')/'repair-check-end-UTC.txt').exists():time.sleep(5)
 c=R/(n+'_R07');normal=(c/'log.checkMesh').read_text();extra=(c/'log.checkMesh-additional').read_text()
 assert '\nEnd\n' in normal and 'Mesh OK.' in normal,(n,'standard quality failure')
 assert '\nEnd\n' in extra and 'Cells with small determinant' not in extra,(n,'determinant failure')
 # Remove only stale copied flag sets when the new native recheck finds zero.
 p=c/'constant/polyMesh/sets/underdeterminedCells'
 if p.exists():p.unlink()
print('ALL THREE MESH QUALITY GATES PASSED',flush=True)
for n in ['B00','C01','C02']:
 print(n,'START20',datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
 p=subprocess.run(['docker','exec','tbm-ramp-1121','bash',str(R/'run_pilot20.sh'),n+'_R07'])
 assert p.returncode==0,(n,'pilot failed; inspect native log')
 p=subprocess.run([sys.executable,str(R/'analyze_pilot20.py'),n]);assert p.returncode==0,(n,'analysis failed; do not erase processor copies')
(R/'PILOTS_COMPLETE.txt').write_text('B00 C01 C02 exactly20 iterations each completed. No 200 iteration run.\n');print('ALL THREE 20 ITERATION PILOTS COMPLETE',flush=True)
