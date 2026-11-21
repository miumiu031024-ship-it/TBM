from pathlib import Path
import shutil,json
R=Path(__file__).resolve().parent
B=Path('/workspace/scratch/tbm-A01-cfd-20261120/B00')
S=R/'SeedSide'
assert not S.exists()
(S/'constant').mkdir(parents=True)
shutil.copytree(B/'constant/polyMesh',S/'constant/polyMesh',ignore=shutil.ignore_patterns('sets'))
for n in ['transportProperties','turbulenceProperties']:shutil.copy2(B/'constant'/n,S/'constant'/n)
shutil.copytree(B/'system',S/'system')
(S/'500').mkdir()
for n in ['U','p','phi','k','omega','nut']:shutil.copy2(B/'500'/n,S/'500'/n)
m=json.loads((B/'case-metadata.json').read_text());m['selected_region']='Water_Xminus';m['parent_final_iteration']=500;m['original_source_case']=str(B)
(S/'case-metadata.json').write_text(json.dumps(m,indent=2))
print('Copied independent seed; original is read-only in container',flush=True)
