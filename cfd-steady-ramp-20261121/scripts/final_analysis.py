"""Postprocess actual ramp experiments. Changes of discretization are not design effects."""
from pathlib import Path
import json,re,csv,hashlib,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-ramp-plots')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-ramp-cache')
Path(os.environ['XDG_CACHE_HOME']).mkdir(parents=True,exist_ok=True)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader,vtkXMLPolyDataReader
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkFiltersCore import vtkProbeFilter,vtkTriangleFilter
from vtkmodules.util.numpy_support import numpy_to_vtk,vtk_to_numpy
from analyze_checkpoint import csv_write
R=Path(__file__).resolve().parent;O=R/'results'

def scalar_patch(path,patch,count):
 b=path.read_bytes();m=re.search(rb'internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(',b);assert m
 st=m.end()+int(m.group(1))*8;st+=re.search(rb'boundaryField\s*\{',b[st:]).end();st+=re.search(rb'\b'+patch.encode()+rb'\s*\{',b[st:]).end()
 m=re.search(rb'value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(',b[st:]);u=re.search(rb'value\s+uniform\s+([\deE+.-]+)\s*;',b[st:])
 if u and (not m or u.start()<m.start()):return np.full(count,float(u.group(1)))
 assert m and m.start()<500 and int(m.group(1))==count
 a=np.frombuffer(b,dtype='<f8',offset=st+m.end(),count=count);assert np.isfinite(a).all();return a

def vtk_final(C):
 n=C.name;series=json.loads((C/'VTK'/f'{n}.vtm.series').read_text())['files'];s=max(series,key=lambda s:s['time']);folder=C/'VTK'/Path(s['name']).stem
 rd=vtkXMLUnstructuredGridReader();rd.SetFileName(str(folder/'internal.vtu'));rd.Update();d=rd.GetOutput();assert d.GetNumberOfCells()==1264925
 return s['time'],folder,rd,d

def sample(data,xyz):
 p=vtkPoints();p.SetData(numpy_to_vtk(xyz,deep=True));pd=vtkPolyData();pd.SetPoints(p)
 q=vtkProbeFilter();q.SetInputData(pd);q.SetSourceData(data);q.Update();out=q.GetOutput()
 mask=vtk_to_numpy(out.GetPointData().GetArray('vtkValidPointMask')).astype(bool);u=vtk_to_numpy(out.GetPointData().GetArray('U')).copy();p=vtk_to_numpy(out.GetPointData().GetArray('p')).copy()*746.200955901462
 assert np.isfinite(u[mask]).all() and np.isfinite(p[mask]).all();return mask,u,p

def process():
 complete=json.loads((R/'SOLVES_COMPLETE.json').read_text());fo=complete['first_order'];so=complete['second_order'];switch=fo['iteration']
 assert fo['normal_solver_end'] and so['normal_solver_end']
 C=R/'F02_second_order_restore';meta=json.loads((C/'case-metadata.json').read_text())
 # Independent provenance validation: identical six actual restart fields.
 hashes={}
 for n in ['U','p','phi','k','omega','nut']:
  a=(R/'F01_upwind_long'/str(switch)/n).read_bytes();b=(C/str(switch)/n).read_bytes();assert a==b,(n,'Restart changed');hashes[n]=hashlib.sha256(a).hexdigest()
 meta['restart_field_sha256']=hashes;meta['restart_hashes_independently_verified']=True;(C/'case-metadata.json').write_text(json.dumps(meta,indent=2))
 fig,axs=plt.subplots(2,3,figsize=(14,7))
 for name,color,label in [('F01_upwind_long','#167ac6','First-order startup'),('F02_second_order_restore','#de6f13','Restored second order')]:
  D=O/name
  for file,col,ax in [('residuals.csv','p',axs[0,0]),('residuals.csv','Ux',axs[0,1]),('boundary-history.csv','outlet_spread_percent',axs[1,0]),('boundary-history.csv','inlet_area_mean_static_Pa',axs[1,1]),('internal-monitor.csv','low_fraction_0p5',axs[1,2])]:
   q=list(csv.DictReader((D/file).open()));x=[float(r['iteration']) for r in q];y=[float(r[col])*(100 if col=='low_fraction_0p5' else .001 if col.endswith('_Pa') else 1) for r in q]
   if file=='residuals.csv':ax.semilogy(x,y,c=color,label=label)
   else:ax.plot(x,y,c=color,label=label)
  q=[json.loads(p.read_text()) for p in D.glob('assessment-*.json')];q.sort(key=lambda x:x['iteration']);axs[0,2].semilogy([s['iteration'] for s in q],[s['native_velocity_L2_last100'] for s in q],'o-',c=color,label=label)
 for ax,label in zip(axs.ravel(),['Initial p residual','Initial Ux residual','Native U relative L2 over last100','6-outlet spread / mean [%]','Area mean inlet static p [kPa]','Monitor fraction speed <0.5 m/s [%]']):
  ax.axvline(switch,c='gray',ls='--',lw=.8);ax.set(xlabel='Steady numerical iteration (not physical time)',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=7)
 axs[0,0].axhline(1e-4,c='black',ls=':',lw=.8)
 axs[0,1].axhline(1e-5,c='black',ls=':',lw=.8)
 axs[0,2].axhline(.01,c='black',ls=':',lw=.8)
 fig.suptitle('Unchanged baseline geometry: long first-order startup then second-order restoration\nNo physical-period or grid-independence inference from this steady history',fontsize=11);fig.tight_layout();fig.savefig(O/'ramp-convergence.png',dpi=180);plt.close(fig)
 cases=[Path('/workspace/scratch/tbm-cause-20261121/S30_upwind'),R/'F01_upwind_long',R/'F02_second_order_restore'];labels=['Short first-order startup','Long first-order startup','Restored second order'];data={};surface={};outlets=[];backflow_history=[];native_ranges={}
 xyz0=np.load('/workspace/scratch/tbm-A01-cfd-20261120/results/B00/surface-r46.4.npz')['points_m'];directions=xyz0/.0464
 for C,label in zip(cases,labels):
  tm,folder,reader,d=vtk_final(C);data[C.name]=(tm,folder,reader,d)
  cu=vtk_to_numpy(d.GetCellData().GetArray('U'));cp=vtk_to_numpy(d.GetCellData().GetArray('p'))*746.200955901462;assert np.isfinite(cu).all() and np.isfinite(cp).all()
  native_ranges[C.name]=dict(maximum_cell_speed_m_s=float(np.linalg.norm(cu,axis=1).max()),minimum_cell_gauge_pressure_Pa=float(cp.min()),maximum_cell_gauge_pressure_Pa=float(cp.max()))
  for radius in [45.2,46.4,47.6]:
   xyz=directions*(radius*.001);xyz=xyz[xyz[:,0]<-.02370];mask,u,p=sample(d,xyz);surface[C.name,radius]=(xyz,mask,u,p)
   np.savez_compressed(O/f'{C.name}-r{radius:.1f}.npz',points_m=xyz,valid=mask,velocity_m_s=u,pressure_Pa=p)
  for i in range(1,7):
   patch=f'outlet{i:02d}';rd=vtkXMLPolyDataReader();rd.SetFileName(str(folder/'boundary'/f'{patch}.vtp'));rd.Update();pd=rd.GetOutput();nf=pd.GetNumberOfCells();assert nf>0
   ids=numpy_to_vtk(np.arange(nf),deep=True);ids.SetName('faceId');pd.GetCellData().AddArray(ids);tri=vtkTriangleFilter();tri.SetInputData(pd);tri.Update();t=tri.GetOutput();v=vtk_to_numpy(t.GetPoints().GetData());f=vtk_to_numpy(t.GetPolys().GetConnectivityArray()).reshape(-1,3);av=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])*.5;assert (av[:,2]<0).all();idx=vtk_to_numpy(t.GetCellData().GetArray('faceId')).astype(int)
   A=np.column_stack([np.bincount(idx,weights=av[:,k],minlength=nf) for k in range(3)]);q=scalar_patch(C/f'{tm:g}'/'phi',patch,nf)*746.200955901462;uv=vtk_to_numpy(pd.GetCellData().GetArray('U'));assert np.max(abs(vtk_to_numpy(pd.GetCellData().GetArray('p'))))<1e-12
   outlets.append(dict(case=C.name,iteration=tm,outlet=patch,native_net_kg_s=float(q.sum()),native_forward_kg_s=float(q[q>0].sum()),native_reverse_kg_s=float(-q[q<0].sum()),minimum_outward_velocity_m_s=float((np.sum(uv*A,axis=1)/np.linalg.norm(A,axis=1)).min()),mass_from_boundary_U_kg_s=float(746.200955901462*np.sum(uv*A))))
   for snap in sorted((p for p in C.iterdir() if p.is_dir() and re.fullmatch(r'\d+(\.\d+)?',p.name) and (p/'phi').exists()),key=lambda p:float(p.name)):
    sq=scalar_patch(snap/'phi',patch,nf)*746.200955901462
    backflow_history.append(dict(case=C.name,iteration=float(snap.name),outlet=patch,native_net_kg_s=float(sq.sum()),native_forward_kg_s=float(sq[sq>0].sum()),native_reverse_kg_s=float(-sq[sq<0].sum())))
 csv_write(O/'outlet-snapshots.csv',outlets)
 csv_write(O/'outlet-backflow-snapshot-history.csv',backflow_history)
 areas=[]
 for radius in [45.2,46.4,47.6]:
  common=np.logical_and.reduce([surface[C.name,radius][1] for C in cases]);ap=2*np.pi*radius**2/100000
  for C in cases:
   xyz,valid,u,p=surface[C.name,radius];speed=np.linalg.norm(u,axis=1)
   for threshold in [.1,.5,1.]:areas.append(dict(case=C.name,iteration=data[C.name][0],radius_mm=radius,threshold_m_s=threshold,common_area_mm2=float(common.sum()*ap),low_area_mm2=float((speed[common]<threshold).sum()*ap),low_fraction=float((speed[common]<threshold).mean()),interpretation='Same geometry, different numerical scheme/startup; not a design improvement'))
 csv_write(O/'low-speed-same-geometry.csv',areas)
 xyz=surface[cases[0].name,46.4][0];common=np.logical_and.reduce([surface[C.name,46.4][1] for C in cases])
 fig=plt.figure(figsize=(13,4.7));axes=[]
 for i,(C,label) in enumerate(zip(cases,labels)):
  _,_,u,p=surface[C.name,46.4];ax=fig.add_subplot(1,3,i+1,projection='3d');axes.append(ax);a=xyz[common]*1000;s=np.linalg.norm(u[common],axis=1);sc=ax.scatter(a[:,0],a[:,1],a[:,2],c=s,cmap='turbo',vmin=0,vmax=1.5,s=1,linewidths=0,depthshade=False)
  ax.set(xlim=(-49,-21),ylim=(-43,43),zlim=(0,43),title=label+f'\niteration {data[C.name][0]:g}');ax.set_box_aspect((28,86,43));ax.view_init(22,-145);ax.set_axis_off()
 fig.colorbar(sc,ax=axes,shrink=.55,label='Speed [m/s]',extend='max');fig.suptitle('Actual CFD on identical interior sphere r=46.4 mm\nSame geometry: changes can reflect numerical diffusion or unsettled fields',fontsize=11);fig.savefig(O/'interior-speed-comparison.png',dpi=180,bbox_inches='tight');plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4))
 for C,label in zip(cases,labels):
  q=[r for r in outlets if r['case']==C.name];ax.plot(range(1,7),[1000*r['native_net_kg_s'] for r in q],'o-',label=label)
 ax.set(xlabel='Left-domain outlet number',ylabel='Native outlet mass flow [g/s]');ax.legend();ax.grid(alpha=.2);fig.suptitle('Saved-state outlet distribution; assess convergence separately');fig.tight_layout();fig.savefig(O/'outlet-comparison.png',dpi=180);plt.close(fig)
 summary=dict(first_order=fo,second_order=so,same_geometry=True,paired_same_restart_verified=True,physical_periodicity_established=False,grid_independence_established=False,low_speed=areas,final_outlets=outlets,native_saved_field_ranges=native_ranges)
 (O/'comparison.json').write_text(json.dumps(summary,indent=2));print('Processed real final fields',flush=True)
if __name__=='__main__':process()
