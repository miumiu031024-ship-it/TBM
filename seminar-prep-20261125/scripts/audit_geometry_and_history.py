"""Read-only source CAD volumes and HISTORICAL flow runtime references."""
from pathlib import Path
import os,csv,re
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-seminar-cache')
import cadquery as cq
repo=Path('/workspace/TBM');out=repo/'seminar-prep-20261125';data=[];base=None
for name in ['B00_baseline','B01_partition_inward_050','C01_partition_yaw_plus1','C02_partition_yaw_minus1']:
 s=cq.importers.importStep(str(repo/'cad-proposals-20261120/water'/f'{name}_water.step')).val()
 parts=sorted(s.Solids(),key=lambda a:a.Center().x);assert len(parts)==3 and s.isValid()
 volumes=[p.Volume(1e-8) for p in parts]
 if base is None:base=volumes
 for reg,p,v,b in zip(['Water_Xminus','Water_Center','Water_Xplus'],parts,volumes,base):
  data.append(dict(model=name[:3],region=reg,CAD_water_volume_mm3=v,volume_change_from_B00_percent=100*(v/b-1),CAD_bounding_zmax_mm=p.BoundingBox().zmax))
with (out/'geometry-volume-by-region.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
rows=[]
for model in ['B00','A01']:
 for filename in ['log.simpleFoam-warmup','log.simpleFoam']:
  p=repo/'a01-cfd-comparison-20261120'/model/filename;t=p.read_text()
  its=[float(v) for v in re.findall(r'^Time = ([\deE+.-]+)',t,re.M)]
  clock=float(re.findall(r'ClockTime = ([\deE+.-]+) s',t)[-1]);assert '\nEnd\n' in t
  rows.append(dict(model=model,log=str(p.relative_to(repo)),first_iteration=min(its),last_iteration=max(its),iterations=len(its),wall_clock_s=clock,wall_seconds_per_iteration=clock/len(its),MPI_processes=4,flow_run_is_historical=True))
with (out/'historical-flow-runtime-reference.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('CAD volumes and historical runtime extracted; no current flow solve')
