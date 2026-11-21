"""Two physical-time pilot calculations; identical water geometry and restart.
A short pilot is not a periodicity or statistically converged mean-flow test.
"""
from pathlib import Path
import shutil,json,re
R=Path(__file__).resolve().parent;seed=R/'SeedSide';base=R/'S00_reference'
for name,dt in [('T30_dt2e5',.00002),('T40_dt1e5',.00001)]:
 C=R/name;assert not C.exists();(C/'constant').mkdir(parents=True);(C/'system').mkdir();(C/'0').mkdir()
 (C/'constant/polyMesh').symlink_to(seed/'constant/polyMesh',target_is_directory=True)
 for n in ['transportProperties','turbulenceProperties']:shutil.copy2(seed/'constant'/n,C/'constant'/n)
 for n in ['U','p','phi','k','omega','nut']:
  b=(base/'500'/n).read_bytes();b=b.replace(b'location    "500";',b'location    "0";');(C/'0'/n).write_bytes(b)
 shutil.copy2(base/'system/decomposeParDict',C/'system/decomposeParDict')
 s=(base/'system/fvSchemes').read_text().replace('default steadyState;','default backward;').replace('bounded Gauss linearUpwind','Gauss linearUpwind').replace('bounded Gauss upwind','Gauss upwind')
 (C/'system/fvSchemes').write_text(s)
 (C/'system/fvSolution').write_text('''FoamFile {version 2.0; format ascii; class dictionary; object fvSolution;}
solvers {
p {solver GAMG; tolerance 1e-8; relTol .01; smoother DICGaussSeidel; cacheAgglomeration true; nCellsInCoarsestLevel 100; directSolveCoarsest yes;}
pFinal {$p; tolerance 1e-7; relTol 0;}
"(U|k|omega)" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol .05;}
"(U|k|omega)Final" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol .01;}
}
PIMPLE {momentumPredictor yes; nOuterCorrectors 2; nCorrectors 1; nNonOrthogonalCorrectors 1;}
relaxationFactors {equations {U .7; k .7; omega .7;}}
''')
 s=(base/'system/controlDict').read_text().replace('application simpleFoam','application pimpleFoam').replace('startFrom latestTime; startTime 500;','startFrom startTime; startTime 0;').replace('endTime 600; deltaT 1;',f'endTime .002; deltaT {dt:.10g}; adjustTimeStep no;')
 s=s.replace('writeControl timeStep; writeInterval 50; purgeWrite 0;','writeControl runTime; writeInterval .001; purgeWrite 3;')
 # Identical physical sample spacing of 0.0001 s for the two time steps.
 interval=round(.0001/dt)
 s=s.replace('executeInterval 5;',f'executeInterval {interval};').replace('writeInterval 5;',f'writeInterval {interval};')
 (C/'system/controlDict').write_text(s)
 (C/'case-metadata.json').write_text(json.dumps(dict(parent_seed='Baseline left-domain steady iteration 500 (unsettled initial field)',cells=1264925,delta_t_s=dt,planned_end_s=.002,time_discretization='backward',outer_correctors=2,pressure_correctors=1,nonorthogonal_correctors=1,physical_sample_interval_s=.0001,purpose='Short URANS pilot and time-step sensitivity; not enough for statistically converged performance or physical periodicity'),indent=2))
print('Prepared physical-time pilot pair; duration 0.002 s each')
