"""Small TBM coolant-domain design variants; no CFD performance claim.

Run: python generate_models.py --source input.step --output OUTPUT
Units: mm. Source STEP is read-only. Requires cadquery2.7/OCP.
"""
import os,json,math,hashlib,argparse,csv
from pathlib import Path
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-cad-cache')
import cadquery as cq
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane,GeomAbs_Sphere,GeomAbs_Cylinder
from OCP.Interface import Interface_Static

def volume(s):return sum(a.Volume(1e-8) for a in s.Solids())
def reg(s):
 b=s.BoundingBox()
 return 'Water_Xplus' if b.xmin>0 else 'Water_Xminus' if b.xmax<0 else 'Water_Center'
def parts(s):return sorted(s.Solids(),key=lambda a:['Water_Xminus','Water_Center','Water_Xplus'].index(reg(a)))
def caps(s):
 out=[]
 for solid in parts(s):
  for f in solid.Faces():
   if f.geomType()=='PLANE' and abs(f.Center().z+10)<1e-5 and abs(f.Area()-math.pi*1.5**2)<1e-4:
    c=f.Center();out.append(dict(region=reg(solid),side='IN' if c.y<0 else 'OUT',x_mm=c.x,y_mm=c.y,z_mm=c.z,area_mm2=f.Area(),radius_mm=1.5))
 for side in ['IN','OUT']:
  for i,p in enumerate(sorted([p for p in out if p['side']==side],key=lambda p:p['x_mm']),1):p['port']=f'{side}_{i:02d}'
 return sorted(out,key=lambda p:p['port'])
def partition_faces(s):
 return [f for f in s.Faces() if f.geomType()=='PLANE' and abs(BRepAdaptor_Surface(f.wrapped).Plane().Axis().Direction().X())>.99999]
def root_edges(s):
 spheres=[f for f in s.Faces() if f.geomType()=='SPHERE'];result=[]
 for f in partition_faces(s):
  for e in f.Edges():
   if any(any(e.isSame(a) for a in q.Edges()) for q in spheres) and not any(e.isSame(a) for a in result):result.append(e)
 return result
def export_water(s,path):
 ass=cq.Assembly(name=path.stem)
 for a in parts(s):ass.add(a,name=reg(a))
 ass.save(str(path),exportType='STEP')
def export_metal(s,path):
 ass=cq.Assembly(name=path.stem)
 ss=sorted(s.Solids(),key=lambda a:a.Center().x)
 for i,a in enumerate(ss):ass.add(a,name='Partition_Xminus' if i==0 else 'Partition_Xplus')
 ass.save(str(path),exportType='STEP')
def read(path):return cq.importers.importStep(str(path)).val()

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--resume',action='store_true');a=parser.parse_args()
 O=a.output
 if O.exists() and not a.resume:raise SystemExit('Refusing to overwrite existing output directory: '+str(O))
 (O/'water').mkdir(parents=True,exist_ok=a.resume);(O/'partition-reference').mkdir(exist_ok=a.resume)
 cache={}
 checkpoint=O/'.validated-models.json'
 if a.resume and checkpoint.exists():
  state=json.loads(checkpoint.read_text());assert state['source_sha256']==hashlib.sha256(a.source.read_bytes()).hexdigest()
  for item in state['models']:
   for p,h in item['hashes'].items():assert hashlib.sha256((O/p).read_bytes()).hexdigest()==h
   cache[item['entry']['id']]=item
 Interface_Static.SetCVal_s('write.step.unit','MM')
 src=read(a.source);ss=src.Solids();fused=ss[0].fuse(*ss[1:]).clean()
 assert fused.isValid() and len(fused.Solids())==3
 baseline=cq.Compound.makeCompound(parts(fused));source_ports=caps(baseline);assert len(source_ports)==36
 positions=sorted({f.Center().x for s in parts(baseline) for f in partition_faces(s)})
 assert len(positions)==4
 centres=[(positions[0]+positions[1])/2,(positions[2]+positions[3])/2]
 assert all(abs((positions[i+1]-positions[i])-1)<1e-6 for i in [0,2])
 shell=cq.Solid.makeSphere(47.9).cut(cq.Solid.makeSphere(44.9)).intersect(cq.Solid.makeBox(120,120,60,cq.Vector(-60,-60,0)))
 fills=[shell.intersect(cq.Solid.makeBox(positions[j+1]-positions[j],120,60,cq.Vector(positions[j],-60,0))) for j in [0,2]]
 pp=parts(baseline);full=pp[0].fuse(*pp[1:],*fills).clean();assert full.isValid() and len(full.Solids())==1
 def slabs_for(shift=0,angle=0):
  slabs=[]
  for c in centres:
   cc=c-math.copysign(shift,c)
   slab=cq.Solid.makeBox(1,140,90,cq.Vector(cc-.5,-70,-20))
   if angle:slab=slab.rotate((cc,0,0),(cc,0,1),angle if c>0 else -angle)
   slabs.append(slab)
  return slabs
 def change(shift=0,angle=0):return full.cut(*slabs_for(shift,angle)).clean()
 rec=change();equivalence=volume(rec.cut(baseline))+volume(baseline.cut(rec));assert equivalence<1e-4
 models=[('B00_baseline','接合済み基準模型',baseline,dict(kind='baseline'))]
 for radius in [.25,.5]:
  s=cq.Compound.makeCompound([x.fillet(radius,root_edges(x)) for x in parts(baseline)])
  assert volume(s.cut(baseline))<1e-5,'Root fillets must add metal, not remove it'
  models.append((f'A0{1 if radius==.25 else 2}_root_R{int(radius*100):03d}',f'隔壁付け根R{radius:g} mm',s,dict(kind='root_fillet',radius_mm=radius)))
 models.append(('B01_partition_inward_050','隔壁2枚を中央側へ各0.5 mm',change(shift=.5),dict(kind='partition_translation',inward_offset_each_mm=.5)))
 for angle in [1,-1]:
  models.append(('C01_partition_yaw_plus1' if angle>0 else 'C02_partition_yaw_minus1',f'隔壁を左右対称に{angle:+d}度傾斜',change(angle=angle),dict(kind='partition_yaw',right_partition_rotation_about_Z_deg=angle,left_partition_rotation_about_Z_deg=-angle,pivot_centres_x_mm=centres)))
 pipes=cq.Compound.makeCompound([cq.Solid.makeCylinder(1.5,10,cq.Vector(p['x_mm'],p['y_mm'],-10),cq.Vector(0,0,1)) for p in source_ports])
 report=[]
 for name,label,s,parameters in models:
  if name in cache:
   item=cache[name];rr=read(O/item['entry']['water_file']);mm=read(O/item['entry']['partition_reference_file'])
   assert rr.isValid() and len(rr.Solids())==3 and mm.isValid() and len(mm.Solids())==2 and len(caps(rr))==36
   report.append(item['entry']);print('Reused hash-verified validated model '+name,flush=True);continue
  assert s.isValid() and len(s.Solids())==3,(name,'invalid water topology')
  current_ports=caps(s);assert len(current_ports)==36
  ports_error=0.
  for p,q in zip(source_ports,current_ports):
   assert p['port']==q['port'] and p['region']==q['region'],(name,p,q)
   err=math.sqrt(sum((p[k]-q[k])**2 for k in ['x_mm','y_mm','z_mm']));ports_error=max(ports_error,err)
   assert err<1e-5 and abs(p['area_mm2']-q['area_mm2'])<1e-5
  order=parts(s);separation=min(order[i].distance(order[j]) for i in range(3) for j in range(i+1,3));assert separation>=1-1e-5
  radial=[];cyl_count=0
  for solid in order:
   for face in solid.Faces():
    ad=BRepAdaptor_Surface(face.wrapped)
    if ad.GetType()==GeomAbs_Sphere and ad.Sphere().Radius()>40:
     sp=ad.Sphere();loc=sp.Location();assert loc.Distance(cq.Vector(0,0,0).toPnt())<1e-6;radial.append(sp.Radius())
    if ad.GetType()==GeomAbs_Cylinder:
     bb=face.BoundingBox();cy=ad.Cylinder()
     if abs(cy.Radius()-1.5)<1e-6 and abs(bb.zmin+10)<1e-5 and abs(bb.zmax)<1e-5:cyl_count+=1
  assert cyl_count==36,(name,cyl_count)
  assert radial and all(min(abs(r-44.9),abs(r-47.9))<1e-6 for r in radial)
  excess=volume(s.cut(full));assert excess<1e-4,(name,'outside unchanged fluid envelope',excess)
  if parameters['kind']=='partition_yaw':
   metal=cq.Compound.makeCompound([shell.intersect(q).clean() for q in slabs_for(angle=parameters['right_partition_rotation_about_Z_deg'])])
  else:metal=full.cut(s).clean()
  assert metal.isValid() and len(metal.Solids())==2,(name,'metal-reference topology',len(metal.Solids()))
  assert volume(metal.intersect(s))<1e-4
  clearance=metal.distance(pipes);assert clearance>.5,(name,'pipe clearance too small',clearance)
  added=volume(s.cut(baseline));removed=volume(baseline.cut(s))
  wf=O/'water'/f'{name}_water.step';mf=O/'partition-reference'/f'{name}_partition_reference.step'
  export_water(s,wf);export_metal(metal,mf)
  imported=read(wf);im=read(mf);assert imported.isValid() and len(imported.Solids())==3 and im.isValid() and len(im.Solids())==2
  mass_delta=abs(volume(imported)-volume(s))
  roundtrip_symmetric_difference=volume(imported.cut(s))+volume(s.cut(imported))
  assert roundtrip_symmetric_difference<1e-4 and mass_delta/volume(s)<1e-5,(name,mass_delta,roundtrip_symmetric_difference)
  imported_ports=caps(imported);assert len(imported_ports)==36
  for p,q in zip(source_ports,imported_ports):assert p['region']==q['region'] and max(abs(p[k]-q[k]) for k in ['x_mm','y_mm','z_mm'])<1e-5
  bb=imported.BoundingBox();bb0=baseline.BoundingBox()
  assert max(abs(getattr(bb,k)-getattr(bb0,k)) for k in ['xmin','xmax','ymin','ymax','zmin','zmax'])<1e-5
  entry=dict(id=name,label=label,parameters=parameters,water_file=str(wf.relative_to(O)),partition_reference_file=str(mf.relative_to(O)),water_solids=3,partition_reference_solids=2,CAD_valid=True,STEP_readback_valid=True,STEP_readback_symmetric_difference_mm3=roundtrip_symmetric_difference,mass_integration_readback_difference_mm3=mass_delta,water_volume_mm3=volume(s),net_water_volume_change_percent=100*(volume(s)/volume(baseline)-1),changed_volume_symmetric_difference_mm3=added+removed,changed_volume_percent=100*(added+removed)/volume(baseline),water_added_mm3=added,water_removed_mm3=removed,minimum_region_separation_mm=separation,minimum_pipe_to_separator_clearance_mm=clearance,port_count=36,inlets=18,outlets=18,maximum_port_position_error_mm=ports_error,unchanged_10mm_long_3mm_bore_cylindrical_faces=cyl_count,outer_radius_mm=47.9,inner_radius_mm=44.9,nominal_radial_water_gap_mm=3,CFD_evaluated=False)
  report.append(entry);print(json.dumps(entry,ensure_ascii=False),flush=True)
  cache[name]=dict(entry=entry,hashes={str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [wf,mf]})
  checkpoint.write_text(json.dumps(dict(source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest(),models=list(cache.values())),indent=2))
 with (O/'ports.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(source_ports[0]));w.writeheader();w.writerows(source_ports)
 data=dict(source_file_name=a.source.name,source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest(),units='mm',representation='Coolant WATER volumes; separate partition-reference files show only divider material within coolant envelope, not complete TBM metal assembly',source_solids=len(ss),baseline_joined_water_regions=3,original_partition_face_x_mm=positions,baseline_reconstruction_symmetric_difference_mm3=equivalence,preserved=['spherical envelope44.9/47.9mm','all36 ports position and3mm bore','straight pipes10mm','minimum divider thickness1mm','three disconnected water domains'],not_evaluated=['CFD velocity/low-speed improvement','thermal/structural performance','manufacturing tolerances, joining method and complete metal assembly'],models=report)
 (O/'geometry-validation.json').write_text(json.dumps(data,indent=2,ensure_ascii=False))
 (O/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(O))+'\n' for p in sorted(O.rglob('*.step'))))
 print('All STEP files exported and reimported successfully',flush=True)
if __name__=='__main__':main()
