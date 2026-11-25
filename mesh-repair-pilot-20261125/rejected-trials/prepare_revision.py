"""New background-grid phase test; source CAD and old meshes stay untouched."""
from pathlib import Path
import sys,shutil,json,re,hashlib,subprocess
R=Path(__file__).resolve().parent
OLD=Path('/workspace/scratch/tbm-seminar-prep-20261125')
REPO=Path('/workspace/TBM')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
if not (R/'protected-before.json').exists():
 files=[REPO/p for p in subprocess.check_output(['git','ls-files','-z'],cwd=REPO).decode().split('\0') if p]
 for model in ['B00','C01','C02']:
  for p in (OLD/model).rglob('*'):
   if p.is_file() and ('polyMesh' in p.parts or p.parent.name in ['system','0'] or p.name.startswith('log.')):files.append(p)
 (R/'protected-before.json').write_text(json.dumps({str(p):sha(p) for p in files},indent=2))
for model in sys.argv[1:] or ['C01']:
 old=OLD/model;C=R/f'{model}_R01';assert not C.exists(),C
 (C/'constant').mkdir(parents=True)
 shutil.copytree(old/'constant/triSurface',C/'constant/triSurface',ignore=shutil.ignore_patterns('*.obj','*.vtk'))
 for f in ['transportProperties','turbulenceProperties']:shutil.copy2(old/'constant'/f,C/'constant'/f)
 shutil.copytree(old/'system',C/'system');shutil.copytree(old/'0',C/'0')
 s=(C/'system/blockMeshDict').read_text()
 m=re.search(r'vertices\s*\((.*?)\);\s*blocks',s,re.S);assert m
 verts=re.findall(r'\(([^()]+)\)',m.group(1));assert len(verts)==8
 shift=.1/28/64/2 # half finest level6 cell width in X; geometry is NOT translated.
 converted=[]
 for v in verts:
  a=list(map(float,v.split()));a[0]+=shift;converted.append('('+' '.join(f'{x:.17g}' for x in a)+')')
 s=s[:m.start(1)]+' '.join(converted)+s[m.end(1):]
 (C/'system/blockMeshDict').write_text(s)
 meta=json.loads((old/'case-metadata.json').read_text())
 meta.update(revision='R01_background_phase_X_half_finest_cell',source_case=str(old),
   background_origin_offset_m=[shift,0,0],CAD_geometry_changed=False,flow_solver_executed=False,
   changed_mesh_parameter='Only the X origin of the background Cartesian grid; sizes/refinement/snap/quality limits unchanged')
 (C/'case-metadata.json').write_text(json.dumps(meta,indent=2))
 (C/'case.foam').touch()
 for p in (C/'constant/triSurface').glob('*'):
  assert sha(p)==sha(old/'constant/triSurface'/p.name),p
 print(model,'R01 ready; grid shift',shift*1000,'mm; original CAD/STL unchanged',flush=True)
