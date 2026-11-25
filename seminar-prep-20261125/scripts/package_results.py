"""Package final seminar/mesh-preparation artifacts only; never run CFD."""
from pathlib import Path
import json,re,shutil,hashlib,zipfile,csv,subprocess,ast
R=Path(__file__).resolve().parent;repo=Path('/workspace/TBM');O=repo/'seminar-prep-20261125'
assert (R/'MESHES_COMPLETE.txt').exists() and (R/'PORT_AUDITS_COMPLETE.txt').exists()
assert 'Memory sampling complete' in (R/'log.memory-monitor').read_text()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
q=json.loads((O/'mesh-quality/mesh-summary.json').read_text());assert len(q)==3
protected=json.loads((R/'protected-before.json').read_text())
for p,h in protected.items():assert sha(Path(p))==h,p
for row in q:
 n=row['model'];C=R/n;D=O/'cases'/n;assert not D.exists();D.mkdir(parents=True)
 assert not (C/'constant/polyMesh').is_symlink() and (C/'constant/polyMesh').resolve().parent==C/'constant'
 assert not list(C.glob('log.simpleFoam*')) and not list(C.glob('log.pimpleFoam*'))
 assert not any(p.is_dir() and p.name!='0' and (p/'U').exists() for p in C.iterdir() if re.fullmatch(r'\d+(?:\.\d+)?',p.name))
 assert 'endTime 0;' in (C/'system/controlDict').read_text()
 assert json.loads((O/f'mesh-quality/{n}-native-region-connectivity.json').read_text())['all36_ports_match_CAD_region']
 shutil.copytree(C/'system',D/'system');shutil.copytree(C/'0',D/'0-initial-conditions')
 (D/'constant-properties').mkdir()
 for f in ['transportProperties','turbulenceProperties']:shutil.copy2(C/'constant'/f,D/'constant-properties'/f)
 for p in C.iterdir():
  if p.is_file() and (p.name.startswith('log.') or p.name in ['case-metadata.json','mesh-start-UTC.txt','mesh-end-UTC.txt','container-cumulative-memory-peak-bytes.txt']):shutil.copy2(p,D/p.name)
 for p in R.glob(f'log.*audit-{n}'):shutil.copy2(p,D/p.name)
S=O/'scripts';S.mkdir(exist_ok=True)
for p in R.glob('*.py'):
 if p.name.startswith('log.'):continue
 ast.parse(p.read_text());shutil.copy2(p,S/p.name)
shutil.copy2(R/'mesh_only.sh',S/'mesh_only.sh')
for f in ['memory-samples.csv','log.prepare','log.mesh-driver','log.port-audit-driver']:
 shutil.copy2(R/f,O/f)
A=O/'audit-logs';A.mkdir(exist_ok=True)
for p in R.glob('log.final-*'):shutil.copy2(p,A/p.name)
for n in ['log.geometry-history-audit','log.high-points','log.outlet-geometry-context','log.finish-preparation','log.slides-final']:
 if (R/n).exists():shutil.copy2(R/n,A/n)
shutil.copy2(repo/'cad-proposals-20261120/ports.csv',O/'CAD-port-map.csv')
P=O/'preservation';P.mkdir(exist_ok=True);shutil.copy2(R/'protected-before.json',P/'pre-existing-file-sha256.json')
(P/'check.json').write_text(json.dumps(dict(pre_existing_files_checked=len(protected),all_protected_hashes_unchanged=True,scope='Pre-existing tracked repository files and the previous review folder; large previous native mesh/field files were not accessed for writing. No claim that every native file was hashed.',CAD_geometry_modified=False,old_mesh_or_result_overwritten=False),indent=2))
(O/'execution-notes.json').write_text(json.dumps(dict(flow_solver_executed=False,first_preflight_failed_before_snappyHexMesh=True,preflight_reason='Container lacked /usr/bin/time; used Bash time keyword instead. No snappy or flow solver started in that first invocation.',all_three_final_mesh_jobs_ended_normally=True,meshes_kept_in_workspace=True,volume_meshes_not_in_download_archive=True,software=dict(OpenFOAM='v2412',image_sha256='5c64f40b9d1f439b931f5c2141ee3b9e1c47f0654722aec4eee1136563652cda',MPI_processes=4,memory_limit_GiB=24),source_repo_commit_before_work=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo).decode().strip()),indent=2))
# Complete validation happens before publication or archive creation.
for p in O.rglob('*.md'):
 for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if not link.startswith(('http:','https:','#')):assert (p.parent/link).exists(),(p,link)
with zipfile.ZipFile(O/'slides/TBM-seminar-2slides.pptx') as z:
 assert z.testzip() is None
 assert len([n for n in z.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml',n)])==2
 assert len([n for n in z.namelist() if re.fullmatch(r'ppt/notesSlides/notesSlide\d+\.xml',n)])==2
from pypdf import PdfReader
pdf=PdfReader(O/'slides/TBM-seminar-2slides.pdf');assert len(pdf.pages)==2
for p in pdf.pages:assert abs(float(p.mediabox.width)/float(p.mediabox.height)-16/9)<1e-5
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
for p in files:
 assert 'polyMesh' not in p.parts and 'triSurface' not in p.parts
 assert p.stat().st_size<40*1024**2,p
 if p.suffix in ['.py','.sh','.md','.json']:
  assert not re.search(rb'(?:ghp_|github_pat_)[A-Za-z0-9_]{20,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----',p.read_bytes()),p
(O/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(O).as_posix()}\n' for p in files))
archive=repo/'TBM-seminar-prep-20261125.zip';assert not archive.exists()
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for p in sorted(O.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(repo))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for n in z.namelist():assert z.read(n)==(repo/n).read_bytes()
print(json.dumps(dict(archive=archive.name,bytes=archive.stat().st_size,sha256=sha(archive),file_count=len(files)+1,protected_files_unchanged=len(protected),flow_solver_executed=False),indent=2))
