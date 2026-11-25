import subprocess,time,csv,datetime,sys
from pathlib import Path
R=Path(__file__).resolve().parent
with (R/'pilot-memory-samples.csv').open('w',buffering=1) as f:
 w=csv.writer(f);w.writerow(['UTC','case','container_memory_current_bytes','container_anon_bytes','container_file_bytes','solver_sum_RSS_KiB','mesh_check_process_count'])
 while not (R/'PILOTS_COMPLETE.txt').exists():
  n=(R/'active-pilot.txt').read_text().strip() if (R/'active-pilot.txt').exists() else 'none'
  p=subprocess.run(['docker','exec','tbm-ramp-1121','bash','-lc',"cat /sys/fs/cgroup/memory.current; awk '$1==\"anon\"||$1==\"file\"{print $1,$2}' /sys/fs/cgroup/memory.stat; ps -eo comm,rss | awk '$1==\"simpleFoam\"{r+=$2}$1==\"checkMesh\"{m++}END{print r+0,m+0}'"],capture_output=True,text=True)
  if p.returncode==0:
   ls=p.stdout.splitlines();w.writerow([datetime.datetime.now(datetime.timezone.utc).isoformat(),n,ls[0],ls[1].split()[1],ls[2].split()[1],*ls[3].split()])
  time.sleep(5)
