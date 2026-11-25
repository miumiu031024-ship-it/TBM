from pathlib import Path
import time,json,subprocess,csv,datetime
R=Path(__file__).resolve().parent
with (R/'memory-samples.csv').open('w',buffering=1) as f:
 w=csv.writer(f);w.writerow(['UTC','case','container_memory_current_bytes','container_memory_peak_bytes'])
 while not (R/'MESHES_COMPLETE.txt').exists() and not (R/'STOP_MONITOR').exists():
  active=next((n for n in ['B00','C01','C02'] if (R/n/'mesh-start-UTC.txt').exists() and not (R/n/'mesh-end-UTC.txt').exists()),'between_cases')
  p=subprocess.run(['docker','exec','tbm-seminar-1125','/bin/bash','-c','cat /sys/fs/cgroup/memory.current /sys/fs/cgroup/memory.peak'],capture_output=True,text=True)
  if p.returncode==0:
   vals=p.stdout.split();w.writerow([datetime.datetime.now(datetime.timezone.utc).isoformat(),active,*vals])
  time.sleep(5)
print('Memory sampling complete')
