from pathlib import Path
import hashlib,json,datetime
R=Path(__file__).resolve().parent;before=json.loads((R/'protected-before.json').read_text());bad=[]
for name,expected in before.items():
 p=Path(name);h=hashlib.sha256()
 if not p.exists():bad.append(dict(file=name,reason='missing'));continue
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 if h.hexdigest()!=expected:bad.append(dict(file=name,reason='hash changed'))
result=dict(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),protected_files_checked=len(before),all_unchanged=not bad,mismatches=bad)
(R/'protected-after-verification.json').write_text(json.dumps(result,indent=2));print(result);assert not bad
