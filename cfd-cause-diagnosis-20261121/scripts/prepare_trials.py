from pathlib import Path
import shutil,json,re,hashlib,numpy as np
R=Path(__file__).resolve().parent; seed=R/'SeedSide'
m=json.loads((seed/'case-metadata.json').read_text()); ports=[p for p in m['ports'] if p['region']=='Water_Xminus']
check=(seed/'log.checkMesh').read_text(); assert 'cells:           1264925' in check or re.search(r'cells:\s+1264925\b',check)
# No cut through an internal face is allowed when isolating a disconnected domain.
boundary=(seed/'constant/polyMesh/boundary').read_text()
ex=re.search(r'oldInternalFaces\s*\{[^}]*nFaces\s+(\d+)',boundary,re.S)
assert ex is None or int(ex.group(1))==0,'Unexpected newly exposed internal surface'
a=np.load('/workspace/scratch/tbm-A01-cfd-20261120/results/B00/surface-r46.4.npz')
xyz=a['points_m'][a['valid'] & (a['points_m'][:,0]<-.02371)][::18]
locations='\n'.join('('+' '.join(f'{x:.14g}' for x in q)+')' for q in xyz)
# Same restart fields and mesh; vary exactly one numerical control per trial.
trials={'S00_reference':{},'S10_pressure1':{'p_relaxation':1},'S20_velocity015':{'U_relaxation':.15},'S30_upwind':{'velocity_convection':'upwind'}}
for name,changes in trials.items():
 c=R/name;assert not c.exists();(c/'constant').mkdir(parents=True)
 (c/'constant/polyMesh').symlink_to(seed/'constant/polyMesh',target_is_directory=True)
 for n in ['transportProperties','turbulenceProperties']:shutil.copy2(seed/'constant'/n,c/'constant'/n)
 (c/'system').mkdir();(c/'500').mkdir()
 for n in ['U','p','phi','k','omega','nut']:shutil.copy2(seed/'500'/n,c/'500'/n)
 # subsetMesh retains zero-face boundary entries. They must not attempt to
 # impose nonzero mass flow on a patch with globally zero area.
 f=c/'500/U';b=f.read_bytes()
 b,count=re.subn(rb'(\binlet(?:0[7-9]|1[0-8])\s*\{\s*type\s+)flowRateInletVelocity',rb'\g<1>zeroGradient',b)
 assert count in (0,12),(c,count)
 f.write_bytes(b)
 for n in ['fvSchemes','fvSolution','decomposeParDict']:shutil.copy2(seed/'system'/n,c/'system'/n)
 f=c/'system/fvSolution';s=f.read_text();s=re.sub(r'residualControl\s*\{[^}]*\}','',s)
 if 'p_relaxation' in changes:s=s.replace('fields {p .2;}','fields {p 1;}')
 if 'U_relaxation' in changes:s=s.replace('equations {U .3;','equations {U .15;')
 f.write_text(s)
 if 'velocity_convection' in changes:
  f=c/'system/fvSchemes';f.write_text(f.read_text().replace('bounded Gauss linearUpwind grad(U)','bounded Gauss upwind'))
 funcs=[]
 for p in ports:
  for suffix,op,field in [('Flux','sum','phi'),('Pressure','areaAverage','p')]:
   n=p['patch'];funcs.append(f'{n}{suffix} {{type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name {n}; operation {op}; fields ({field}); writeFields false; log false; writeControl timeStep; writeInterval 1;}}')
 funcs.append('vort {type vorticity; libs (fieldFunctionObjects); result vorticity; executeControl timeStep; executeInterval 5; writeControl none;}')
 funcs.append('interiorProbes {type probes; libs (sampling); fields (U p k omega vorticity); interpolationScheme cellPoint; executeControl timeStep; executeInterval 5; writeControl timeStep; writeInterval 5; probeLocations (\n'+locations+'\n);}')
 control='''FoamFile {version 2.0; format ascii; class dictionary; object controlDict;}
application simpleFoam; startFrom latestTime; startTime 500; stopAt endTime; endTime 600; deltaT 1;
writeControl timeStep; writeInterval 50; purgeWrite 0; writeFormat binary; writePrecision 10; writeCompression off;
timeFormat general; timePrecision 10; runTimeModifiable true;
functions {\n'''+ '\n'.join(funcs)+'\n}\n'
 (c/'system/controlDict').write_text(control)
 (c/'case-metadata.json').write_text(json.dumps(dict(parent=m,changes=changes,original_region_cell_count=1264925,new_internal_boundary_faces=0,initial_iteration=500,end_iteration=600,monitor_points=len(xyz),purpose='One-factor numerical sensitivity, not physical time'),indent=2))
print('Prepared',len(trials),'one-factor cases, each identical 1,264,925-cell side domain; probes:',len(xyz))
(R/'trials.json').write_text(json.dumps(trials,indent=2))
