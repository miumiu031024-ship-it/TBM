"""Prepare NEW B00/C01/C02 mesh-only cases from existing immutable STEP files."""
from pathlib import Path
import hashlib,json,subprocess
R=Path(__file__).resolve().parent
repo=Path('/workspace/TBM')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
tracked=subprocess.check_output(['git','ls-files','-z'],cwd=repo).decode().split('\0')
protected={str(repo/n):sha(repo/n) for n in tracked if n and (repo/n).is_file()}
for p in (repo/'structure-review-20261124').rglob('*'):
 if p.is_file():protected[str(p)]=sha(p)
(R/'protected-before.json').write_text(json.dumps(protected,indent=2))
source=(repo/'a01-cfd-comparison-20261120/prepare.py').read_text()
source=source.replace("for name,file in [('B00','B00_baseline_water.step'),('A01','A01_root_R025_water.step')]:", "for name,file in [('B00','B00_baseline_water.step'),('C01','C01_partition_yaw_plus1_water.step'),('C02','C02_partition_yaw_minus1_water.step')]:")
# Generalize the root selector to the existing tilted planar dividers, in memory only.
source=source.replace('g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)', '''g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
from OCP.BRepAdaptor import BRepAdaptor_Surface
g.partition_faces=lambda s:[f for f in s.Faces() if f.geomType()=='PLANE' and abs(BRepAdaptor_Surface(f.wrapped).Plane().Axis().Direction().X())>.95]''')
source=source.replace("(C/'constant/triSurface/rootGuide.stl').write_text(guide)", '''# Same distance-refinement rule around each model's actual divider roots.
    localguide=['solid rootGuide\\n'];edge_count=0
    for solid in g.parts(s):
        for e in g.root_edges(solid):
            edge_count+=1;n=max(2,math.ceil(e.Length()/.10))
            vv=np.array([e.positionAt(t).toTuple() for t in np.linspace(0,1,n)])*.001
            d=np.array([.000001,0,0])
            for aa,bb in zip(vv[:-1],vv[1:]):
                localguide.extend([facet(np.array([aa-d,aa+d,bb+d])),facet(np.array([aa-d,bb+d,bb-d]))])
    assert edge_count>=8,(name,edge_count)
    localguide.append('endsolid rootGuide\\n')
    (C/'constant/triSurface/rootGuide.stl').write_text(''.join(localguide))''')
# Common monitoring sites omit the divider slabs of all three comparison shapes.
old='for lo,hi in zip(planes[::2],planes[1::2]):keep&=~((mx>lo-.01)&(mx<hi+.01))'
new='''for lo,hi in zip(planes[::2],planes[1::2]):keep&=~((mx>lo-.01)&(mx<hi+.01))
for angle in [1,-1]:
    for centre in [-23.199977677315708,23.199999069163503]:
        a=np.deg2rad(angle if centre>0 else -angle)
        distance=(mx-centre)*np.cos(a)+monitor_xyz[:,1]*1000*np.sin(a)
        keep &= np.abs(distance)>.51'''
assert old in source;source=source.replace(old,new)
source=source.replace('endTime 500; deltaT 1;', 'endTime 0; deltaT 1;')
source=source.replace("mf.write_text(json.dumps(md,indent=2))",'''md['stage']='mesh_preparation_only_no_flow_solve_authorized_in_this_job'
    md['root_refinement']='same distances and levels applied to actual roots of each existing model'
    md['root_guide_edge_count']=edge_count
    md['flow_solver_executed']=False
    md['proposed_pilot']='200 iterations per model:100 upwind then100 linearUpwind; not executed'
    mf.write_text(json.dumps(md,indent=2))''')
compile(source,str(R/'prepare_new.py'),'exec')
(R/'prepare_new.py').write_text(source)
print('Prepared new generator; existing repository files preserved',flush=True)
