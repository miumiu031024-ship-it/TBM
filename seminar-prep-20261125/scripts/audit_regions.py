"""Check native face owners against the three named cellZones. Read-only."""
from pathlib import Path
import re,sys,json
import numpy as np
R=Path(__file__).resolve().parent;D=Path('/workspace/TBM/seminar-prep-20261125/mesh-quality');D.mkdir(exist_ok=True)
names=['Water_Xminus','Water_Center','Water_Xplus']
for n in sys.argv[1:] or ['B00','C01','C02']:
 C=R/n;assert (C/'mesh-end-UTC.txt').exists()
 mesh=C/'constant/polyMesh';b=(mesh/'owner').read_bytes();assert b'LSB;label=32' in b[:1024]
 nc=int(re.search(rb'nCells:(\d+)',b[:1200]).group(1));m=re.search(rb'\n(\d+)\n\(',b);assert m
 owner=np.frombuffer(b,dtype='<i4',count=int(m.group(1)),offset=m.end());assert owner.min()>=0 and owner.max()<nc
 z=(mesh/'cellZones').read_bytes();membership=np.full(nc,-1,dtype=np.int8);counts={}
 for name in names:
  m=re.search(rb'\b'+name.encode()+rb'\s*\{\s*type\s+cellZone;\s*cellLabels\s+List<label>\s+(\d+)\s*\(',z);assert m,name
  count=int(m.group(1));start=m.end();labels=np.frombuffer(z,dtype='<i4',count=count,offset=start)
  assert labels.min()>=0 and labels.max()<nc and (membership[labels]<0).all()
  membership[labels]=names.index(name);counts[name]=count
 assert (membership>=0).all() and sum(counts.values())==nc
 boundary=(mesh/'boundary').read_text();info={a:(int(c),int(s)) for a,c,s in re.findall(r'(\w+)\s*\{[^{}]*?nFaces\s+(\d+);[^{}]*?startFace\s+(\d+);[^{}]*?\}',boundary)}
 meta=json.loads((C/'case-metadata.json').read_text());ports=[]
 for p in meta['ports']:
  count,start=info[p['patch']];belong=membership[owner[start:start+count]]
  assert count>0 and (belong==names.index(p['region'])).all(),p['patch']
  ports.append(dict(patch=p['patch'],CAD_port=p['port'],region=p['region'],native_faces=count,all_faces_match_expected_region=True))
 for region in names:
  pp=[p for p in ports if p['region']==region];assert len(pp)==12
  assert sum(p['patch'].startswith('inlet') for p in pp)==6
 result=dict(model=n,cells=nc,cell_zone_counts=counts,all_cells_in_exactly_one_named_zone=True,all36_ports_match_CAD_region=True,ports=ports,flow_solver_executed=False)
 (D/f'{n}-native-region-connectivity.json').write_text(json.dumps(result,indent=2))
 print(n,counts,'all36 ports match intended domains',flush=True)
