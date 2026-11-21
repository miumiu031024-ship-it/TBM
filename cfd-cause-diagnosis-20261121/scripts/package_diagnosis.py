"""Validate completed investigations, then package and copy into repo. No git writes."""
from pathlib import Path
import json,re,shutil,csv,hashlib,zipfile,sys,platform
import numpy as np
R=Path(__file__).resolve().parent;O=R/'diagnosis';repo=Path('/workspace/TBM')
steady=json.loads((O/'steady-trials-summary.json').read_text());transient=json.loads((O/'transient-pilots-summary.json').read_text())
assert len(steady)==4
for n,s in steady.items():assert s['completed_normally'] and s['last_iteration']==600 and not s['residual_goals_met'],n
for n,dt in [('T30_dt2e5',2e-5),('T40_dt1e5',1e-5)]:
 s=transient[n];assert s['completed_normally'] and abs(s['last_physical_time_s']-.002)<1e-10 and s['delta_t_s']==dt
 assert s['completed_steps']==round(.002/dt)
 assert s['maximum_courant_number']<2
 assert abs(s['final_flows']['mass_imbalance_fraction'])<1e-4
 d=np.load(O/f'{n}-physical-probes.npz');assert d['velocity_m_s'].shape==(20,1359,3)
 assert np.allclose(d['time_s'],np.arange(1,21)*.0001)
 assert d['valid'].all()
 assert np.isfinite(d['velocity_m_s']).all() and np.isfinite(d['pressure_Pa']).all() and np.isfinite(d['vorticity_s_inv']).all()
assert 'paired_final_velocity_relative_L2' in transient
# Exact same restart U in all one-factor tests, exact same U in the transient pair.
sh=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for field in ['U','p','phi','k','omega','nut']:
 assert len({sh(R/n/'500'/field) for n in steady})==1
 assert len({sh(R/n/'0'/field) for n in ['T30_dt2e5','T40_dt1e5']})==1
seed=R/'SeedSide';mesh=seed/'constant/polyMesh'
names=[*steady,'T30_dt2e5','T40_dt1e5']
assert all((R/n/'constant/polyMesh').resolve()==mesh.resolve() for n in names)
for property_name in ['transportProperties','turbulenceProperties']:
 assert len({sh(R/n/'constant'/property_name) for n in names})==1
boundary=(mesh/'boundary').read_text()
exposed=re.search(r'oldInternalFaces\s*\{[^}]*nFaces\s+(\d+)',boundary,re.S)
assert exposed is None or int(exposed.group(1))==0
checks=dict(steady_tests_completed=4,physical_pilots_completed=2,physical_duration_s=.002,physical_probe_frames_each=20,identical_steady_restart_U=True,identical_physical_restart_U=True,new_internal_faces_from_subsetting=0,selected_domain_cells=1264925,periodicity_established=False,statistical_mean_convergence_established=False,A01_improvement_evaluated_in_this_new_study=False)
checks['shared_mesh_file_hashes']={p.name:sh(p) for p in mesh.iterdir() if p.is_file()}
checks['original_baseline_metadata']=json.loads((seed/'case-metadata.json').read_text())
(O/'validation-and-provenance.json').write_text(json.dumps(checks,indent=2))
# Save exact controls, metadata, raw monitor histories, and solver logs.
for name in [*steady,'T30_dt2e5','T40_dt1e5']:
 C=R/name;D=O/'cases'/name;D.mkdir(parents=True,exist_ok=True)
 shutil.copytree(C/'system',D/'system',dirs_exist_ok=True)
 (D/'properties').mkdir(exist_ok=True)
 for n in ['transportProperties','turbulenceProperties']:shutil.copy2(C/'constant'/n,D/'properties'/n)
 shutil.copy2(C/'case-metadata.json',D/'case-metadata.json')
 shutil.copytree(C/'postProcessing',D/'raw-postProcessing',dirs_exist_ok=True)
 for p in C.glob('log.*'):shutil.copy2(p,D/p.name)
D=O/'seed-metadata';D.mkdir(exist_ok=True)
for n in ['log.subsetMesh','log.checkMesh','case-metadata.json']:shutil.copy2(seed/n,D/n)
D=O/'baseline-BC-reference';D.mkdir(exist_ok=True)
for n in ['U','p','k','omega','nut']:shutil.copy2(Path('/workspace/scratch/tbm-A01-cfd-20261120/B00/0')/n,D/n)
(D/'README.txt').write_text('Original baseline boundary-condition reference at initial zero field. Actual tests restart from baseline iteration 500, whose full volume fields are excluded from the download. Zero-face inactive inlet entries are changed to zeroGradient after subsetting; all six active inlet mass flow rates are unchanged.\n')
D=O/'scripts';D.mkdir(exist_ok=True)
for p in sorted(R.glob('*.py')):
 if p.name!='queue_transient.py':shutil.copy2(p,D/p.name)
for n in ['run_trials.sh','run_transient.sh']:shutil.copy2(R/n,D/n)
(O/'software-context.json').write_text(json.dumps(dict(Python=sys.version,platform=platform.platform(),OpenFOAM='v2412',image_id='sha256:5c64f40b9d1f439b931f5c2141ee3b9e1c47f0654722aec4eee1136563652cda',MPI_processes=4,CPU_quota=4),indent=2))
(O/'preflight-not-used.json').write_text(json.dumps(dict(case='T10_dt1e4',delta_t_s=1e-4,intentionally_stopped=True,reason='Maximum Courant number approximately 5 in startup; refined time steps for the final pilot pair. Not a solver blow-up or a final study result.',included_in_comparison=False),indent=2))
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
for p in files:
 assert p.stat().st_size<30*1024**2,(p,p.stat().st_size)
 assert not any(part in p.parts for part in ['polyMesh','processor0','preflight-high-Courant'])
(O/'SHA256SUMS.txt').write_text(''.join(f'{sh(p)}  {p.relative_to(O).as_posix()}\n' for p in files))
dest=repo/'cfd-cause-diagnosis-20261121';archive=repo/'TBM-CFD-cause-diagnosis-20261121.zip'
assert not dest.exists() and not archive.exists()
shutil.copytree(O,dest)
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for p in sorted(dest.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(repo))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for n in z.namelist():assert z.read(n)==(repo/n).read_bytes()
root=repo/'README.md';s=root.read_text();first,tail=s.split('\n',1)
intro='''
**2026-11-21：未収束・出口流量偏りの原因切り分けを追加しました。** 出口管内部の局所的な逆向き流れと数値設定への感度を確認しました。周期性・単一の原因・A01の改善はまだ確定していません。

- [原因切り分けの図・追加計算・数値と限界](cfd-cause-diagnosis-20261121/README.md)
- [原因調査データZIPをダウンロード](https://github.com/miumiu031024-ship-it/TBM/raw/refs/heads/main/TBM-CFD-cause-diagnosis-20261121.zip)

'''
root.write_text(first+'\n'+intro+tail.lstrip('\n'))
print(json.dumps(dict(archive=archive.name,bytes=archive.stat().st_size,sha256=sh(archive),packaged_files=len(files)+1),indent=2))
