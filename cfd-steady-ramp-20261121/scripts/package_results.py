"""Package reviewed ramp outputs; preserve previous archives and unrelated files."""
from pathlib import Path
import json,re,hashlib,shutil,zipfile,csv,sys,platform
import numpy as np
R=Path(__file__).resolve().parent;O=R/'results';repo=Path('/workspace/TBM');s=json.loads((O/'comparison.json').read_text());fo=s['first_order'];so=s['second_order'];plan=json.loads((R/'plan.json').read_text())
assert fo['normal_solver_end'] and so['normal_solver_end']
assert 1000<=fo['iteration']<=2000 and 400<=so['iteration']-fo['iteration']<=600
assert s['paired_same_restart_verified'] and s['same_geometry'] and not s['physical_periodicity_established']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name,q in [('F01_upwind_long',fo),('F02_second_order_restore',so)]:
 C=R/name;M=json.loads((C/'case-metadata.json').read_text());assert M['cells']==1264925
 assert (C/'constant/polyMesh').resolve()==Path('/workspace/scratch/tbm-cause-20261121/SeedSide/constant/polyMesh').resolve()
 logs=list(C.glob('log.simpleFoam-*'));assert logs and all('\nEnd\n' in p.read_text() for p in logs)
 for radius in [45.2,46.4,47.6]:
  d=np.load(O/f'{name}-r{radius:.1f}.npz');mask=d['valid'];assert mask.sum()>20000
  assert np.isfinite(d['velocity_m_s'][mask]).all() and np.isfinite(d['pressure_Pa'][mask]).all()
 assert q['all_steady_numerical_criteria_met']==all(q['criteria'].values())
 D=O/'cases'/name;D.mkdir(parents=True,exist_ok=True);shutil.copytree(C/'system',D/'system',dirs_exist_ok=True)
 (D/'properties').mkdir(exist_ok=True)
 for n in ['transportProperties','turbulenceProperties']:shutil.copy2(C/'constant'/n,D/'properties'/n)
 shutil.copy2(C/'case-metadata.json',D/'case-metadata.json');shutil.copytree(C/'postProcessing',D/'raw-postProcessing',dirs_exist_ok=True)
 for p in C.glob('log.*'):shutil.copy2(p,D/p.name)
for rank in range(4):
 a=R/f'F01_upwind_long/processor{rank}/constant/polyMesh/cellProcAddressing';b=R/f'F02_second_order_restore/processor{rank}/constant/polyMesh/cellProcAddressing';assert sha(a)==sha(b)
for field in ['U','p','phi','k','omega','nut']:
 a=R/'F01_upwind_long'/str(fo['iteration'])/field;b=R/'F02_second_order_restore'/str(fo['iteration'])/field;assert sha(a)==sha(b)
 assert json.loads((R/'F02_second_order_restore/case-metadata.json').read_text())['restart_field_sha256'][field]==sha(b)
for prop in ['transportProperties','turbulenceProperties']:assert sha(R/'F01_upwind_long/constant'/prop)==sha(R/'F02_second_order_restore/constant'/prop)
assert (R/'F01_upwind_long/system/fvSolution').read_bytes()==(R/'F02_second_order_restore/system/fvSolution').read_bytes()
a=(R/'F01_upwind_long/system/fvSchemes').read_text();b=(R/'F02_second_order_restore/system/fvSchemes').read_text();assert a.replace('div(phi,U) bounded Gauss upwind;','div(phi,U) bounded Gauss linearUpwind grad(U);')==b
loc=list(csv.DictReader((O/'spatial-update-location.csv').open()))
for name,q in [('F01_upwind_long',fo),('F02_second_order_restore',so)]:assert any(r['case']==name and int(r['iteration'])==q['iteration'] for r in loc)
outs=list(csv.DictReader((O/'outlet-snapshots.csv').open()))
for name,q in [('F01_upwind_long',fo),('F02_second_order_restore',so)]:
 rows=[r for r in outs if r['case']==name];assert len(rows)==6
 assert abs(sum(float(r['native_net_kg_s']) for r in rows)-q['final_monitors']['outlet_kg_s'])<1e-8
for n in ['PLAN_JA.md','plan.json']:shutil.copy2(R/n,O/n)
D=O/'scripts';D.mkdir(exist_ok=True)
for p in R.glob('*.py'):shutil.copy2(p,D/p.name)
shutil.copy2(R/'run_block.sh',D/'run_block.sh')
# Reference initial BC definitions; the actual full restart fields stay in the workspace.
D=O/'baseline-BC-reference';D.mkdir(exist_ok=True)
for n in ['U','p','k','omega','nut']:shutil.copy2(Path('/workspace/scratch/tbm-A01-cfd-20261120/B00/0')/n,D/n)
(D/'README.txt').write_text('Original baseline initial BC reference, not the full600/phase-switch restart fields. In the selected left domain only six inlets and six outlets have faces. Zero-face unused inlet patches are inert zeroGradient. Native restart hashes and actual controls are provided in case metadata.\n')
(O/'execution-notes.json').write_text(json.dumps(dict(reconstruction_after_1400_recovered=True,solver_1200_to_1400_ended_normally=True,solver_steps_rerun_due_to_recovery=False,reason='On-disk shell runner was edited while waiting; resumed shell offsets failed only after normal solverEnd. Reconstructed saved1300/1400 fields and continued. Later runners are immutable in-memory scripts.'),indent=2))
(O/'validation.json').write_text(json.dumps(dict(normally_ended_blocks=sum(len(list((R/n).glob('log.simpleFoam-*'))) for n in ['F01_upwind_long','F02_second_order_restore']),cells=1264925,same_geometry_mesh_properties_BCs=True,restart_six_fields_equal=True,native_face_flux_matches_monitors=True,first_order_numeric_criteria=fo['all_steady_numerical_criteria_met'],second_order_numeric_criteria=so['all_steady_numerical_criteria_met'],physical_period_established=False,grid_independence_established=False,A01_design_effect_evaluated=False,Python=sys.version,platform=platform.platform(),OpenFOAM='v2412',MPI_processes=4),indent=2))
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
for p in files:assert p.stat().st_size<50*1024**2 and 'polyMesh' not in p.parts and not p.name.endswith('.step'),p
(O/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(O).as_posix()}\n' for p in files))
dest=repo/'cfd-steady-ramp-20261121';archive=repo/'TBM-CFD-steady-ramp-20261121.zip';assert not dest.exists() and not archive.exists()
shutil.copytree(O,dest)
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for p in sorted(dest.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(repo))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for n in z.namelist():assert z.read(n)==(repo/n).read_bytes()
root=repo/'README.md';body=root.read_text();first,tail=body.split('\n',1)
verdict='二次精度でも今回の数値収束基準を満たしました。格子独立性は未検証です。' if so['all_steady_numerical_criteria_met'] else '今回の上限までに二次精度での数値収束基準を満たせませんでした。'
intro=f'''
**2026-11-21：一次精度の長い立上げ→二次精度への復帰を実行しました。** {verdict}

- [立上げ・復帰の実測結果と図](cfd-steady-ramp-20261121/README.md)
- [計算履歴・条件・結果ZIPをダウンロード](https://github.com/miumiu031024-ship-it/TBM/raw/refs/heads/main/TBM-CFD-steady-ramp-20261121.zip)

'''
root.write_text(first+'\n'+intro+tail.lstrip('\n'))
print(json.dumps(dict(archive=archive.name,bytes=archive.stat().st_size,sha256=sha(archive),file_count=len(files)+1),indent=2))
