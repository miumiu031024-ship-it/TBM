"""Checkpointed first-order warmup, followed by second-order restoration.
No geometry changes and no physical-time interpretation of steady iterations.
"""
from pathlib import Path
import subprocess,time,json,shutil,os,hashlib
from analyze_checkpoint import assess
R=Path(__file__).resolve().parent;P=json.loads((R/'plan.json').read_text())

def state(**kw):(R/'driver-state.json').write_text(json.dumps(kw,indent=2))
def block(name,start,end):
 state(phase=name,start=start,target=end,status='solving_or_recovering')
 C=R/name;log=C/f'log.simpleFoam-{start}-{end}'
 if log.exists():
  assert '\nEnd\n' in log.read_text(),f'Incomplete solver log requires explicit recovery: {log}'
  rec=C/f'log.reconstruct-{end}'
  if not rec.exists() or '\nEnd\n' not in rec.read_text():
   script=f'source /usr/lib/openfoam/openfoam2412/etc/bashrc; cd {C}; reconstructPar -time "{end-100},{end}" -fields "(U p phi k omega nut)" > log.reconstruct-{end} 2>&1'
   subprocess.run(['docker','exec','tbm-ramp-1121','bash','-lc',script],check=True)
 else:
  # Freeze the complete runner as an in-memory bash -c argument. Editing the
  # source file later cannot change a running shell's file offsets.
  script=(R/'run_block.sh').read_text();(C/f'log.block-script-{end}.sh').write_text(script)
  with (R/f'log.driver-{name}-{end}').open('w') as f:
   subprocess.run(['docker','exec','tbm-ramp-1121','bash','-c',script,'ramp-block',name,str(start),str(end)],cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
 return assess(name,end)

state(phase='F01_upwind_long',start=600,target=800,status='waiting_initial_block')
for i in range(1800):
 p=R/'F01_upwind_long/log.reconstruct-800'
 if p.exists() and '\nEnd\n' in p.read_text():break
 q=R/'F01_upwind_long/log.simpleFoam-600-800'
 if q.exists() and ('FOAM FATAL' in q.read_text() or 'Floating point exception (8)' in q.read_text()):raise RuntimeError('Initial solver block failed')
 time.sleep(2)
else:raise TimeoutError('Initial block did not complete')
summary=assess('F01_upwind_long',800);end=800
while end<P['first_order']['maximum_iteration']:
 if end>=1000 and summary['all_steady_numerical_criteria_met']:break
 summary=block('F01_upwind_long',end,end+200);end+=200
fo_summary=summary
state(phase='first_order_complete',final=end,criteria_met=summary['all_steady_numerical_criteria_met'])
# Restore only the velocity convection scheme. All other numerical controls,
# mesh, BCs, turbulence model, properties and restart field remain unchanged.
S=R/'F01_upwind_long';C=R/'F02_second_order_restore';assert not C.exists()
(C/'constant').mkdir(parents=True);(C/'system').mkdir();(C/str(end)).mkdir()
(C/'constant/polyMesh').symlink_to((S/'constant/polyMesh').resolve(),target_is_directory=True)
for n in ['transportProperties','turbulenceProperties']:shutil.copy2(S/'constant'/n,C/'constant'/n)
for n in ['U','p','phi','k','omega','nut']:shutil.copy2(S/str(end)/n,C/str(end)/n)
for n in ['fvSchemes','fvSolution','controlDict','decomposeParDict']:shutil.copy2(S/'system'/n,C/'system'/n)
p=C/'system/fvSchemes';s=p.read_text();assert s.count('div(phi,U) bounded Gauss upwind;')==1;p.write_text(s.replace('div(phi,U) bounded Gauss upwind;','div(phi,U) bounded Gauss linearUpwind grad(U);'))
p=C/'system/controlDict';s=p.read_text().replace('startTime 600;',f'startTime {end};');p.write_text(s)
meta=json.loads((S/'case-metadata.json').read_text());meta['restart_field_sha256']={n:hashlib.sha256((C/str(end)/n).read_bytes()).hexdigest() for n in ['U','p','phi','k','omega','nut']};meta.update(parent=str(S),restart_iteration=end,first_order_criteria_met_at_restart=summary['all_steady_numerical_criteria_met'],changed_setting='velocity_convection_upwind_to_linearUpwind_only');(C/'case-metadata.json').write_text(json.dumps(meta,indent=2))
start=end
for target in range(start+200,start+P['second_order']['maximum_extra_iterations']+1,200):
 summary=block(C.name,end,target);end=target
 if end>=start+400 and summary['all_steady_numerical_criteria_met']:break
state(phase='complete',first_order=fo_summary,second_order=summary)
# Native final fields are already reconstructed; convert actual saved states.
for name in ['F01_upwind_long','F02_second_order_restore']:
 script=f'source /usr/lib/openfoam/openfoam2412/etc/bashrc; cd /workspace/scratch/tbm-steady-ramp-20261121/{name}; foamToVTK -latestTime -fields "(U p)" > log.foamToVTK 2>&1'
 subprocess.run(['docker','exec','tbm-ramp-1121','bash','-lc',script],check=True)
(R/'SOLVES_COMPLETE.json').write_text(json.dumps(dict(first_order=fo_summary,second_order=summary),indent=2))
print('Both stages completed and real final fields exported',flush=True)
