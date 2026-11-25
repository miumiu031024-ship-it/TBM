from pathlib import Path
import re,json,sys,hashlib
R=Path(__file__).resolve().parent
for n in sys.argv[1:]:
 c=R/(n+'_R07');assert (c/'case-metadata.json').exists();assert not (c/'log.simpleFoam20').exists()
 p=c/'system/controlDict';s=p.read_text().replace('endTime 0;','endTime 20;').replace('writeInterval 100;','writeInterval 20;').replace('writeCompression off;','writeCompression on;').replace('writeInterval 10;','writeInterval 1;').replace('executeInterval 20;','executeInterval 1;')
 # Probe output every iteration, whole volume only at 20.
 s=s.replace('executeInterval 1; writeControl timeStep; writeInterval 20;','executeInterval 1; writeControl timeStep; writeInterval 1;')
 s=s.replace('runTimeModifiable true;','runTimeModifiable false;')
 p.write_text(s)
 p=c/'system/fvSolution';s=p.read_text();s=re.sub(r'residualControl\s*\{[^{}]*\}','',s);p.write_text(s)
 p=c/'system/fvSchemes';s=p.read_text().replace('div(phi,U) bounded Gauss linearUpwind grad(U);','div(phi,U) bounded Gauss upwind;');assert 'div(phi,U) bounded Gauss upwind;' in s;p.write_text(s)
 meta=json.loads((c/'case-metadata.json').read_text());meta.update(pilot=dict(iterations=20,initialization='identical original uniform initial fields; no potentialFoam initialization',convection_U='bounded Gauss upwind for startup only',startTime=0,endTime=20,compress_output=True,parallel_ranks=4,functions='36 patch flux and mean pressure each iteration, 7864 internal U probes each iteration, yPlus at writeTime',purpose='startup stability and runtime/memory measurement, NOT design improvement/convergence validation'))
 (c/'case-metadata.json').write_text(json.dumps(meta,indent=2));print(n,'configured EXACT 20 upwind cold-start iterations',flush=True)
