"""Audit real matched CFD runs. No analytic/illustrative substitute for CFD."""
from pathlib import Path
import os,json,re,csv,math,sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-a01-plots')
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader,vtkXMLPolyDataReader
from vtkmodules.vtkFiltersCore import vtkProbeFilter,vtkTriangleFilter
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.util.numpy_support import vtk_to_numpy,numpy_to_vtk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from probe_history import analyze_probe_history

R=Path(__file__).resolve().parent
O=R/'results';O.mkdir(exist_ok=True)

def scalar_patch(path,patch,count):
    b=path.read_bytes()
    assert b'LSB' in b[:1024] and b'scalar=64' in b[:1024]
    m=re.search(rb'internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(',b)
    start=m.end()+int(m.group(1))*8 if m else 0
    start+=re.search(rb'boundaryField\s*\{',b[start:]).end()
    m=re.search(rb'\b'+patch.encode()+rb'\s*\{',b[start:]);assert m,patch;start+=m.end()
    m=re.search(rb'value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(',b[start:])
    u=re.search(rb'value\s+uniform\s+([\deE+.-]+)\s*;',b[start:])
    if u and (not m or u.start()<m.start()):return np.full(count,float(u.group(1)))
    assert m and m.start()<500 and int(m.group(1))==count,(path,patch,count)
    a=np.frombuffer(b,dtype='<f8',offset=start+m.end(),count=count);assert np.isfinite(a).all();return a

def native_u(path):
    b=path.read_bytes();m=re.search(rb'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(',b);assert m
    a=np.frombuffer(b,dtype='<f8',offset=m.end(),count=3*int(m.group(1))).reshape(-1,3)
    assert np.isfinite(a).all();return a

def csv_write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def monitors(C):
    allrows={};port_areas={}
    for port in json.loads((C/'case-metadata.json').read_text())['ports']:
        for suffix in ['Flux','Pressure']:
            rows={}
            for p in sorted((C/'postProcessing'/(port['patch']+suffix)).glob('*/surfaceFieldValue.dat'),key=lambda p:float(p.parent.name)):
                body=p.read_text();am=re.search(r'#\s*Area\s*:\s*([\deE+.-]+)',body)
                assert am,(p,'Missing native patch area')
                area=float(am.group(1));assert area>0
                if port['patch'] in port_areas:assert abs(port_areas[port['patch']]-area)<1e-12
                port_areas[port['patch']]=area
                for line in body.splitlines():
                    if line.strip() and not line.lstrip().startswith('#'):
                        a=list(map(float,line.split()));rows[a[0]]=a[1]
            assert rows,(C,port,suffix)
            allrows[port['patch']+suffix]=rows
    return allrows,port_areas

def probe(data,xyz):
    pts=vtkPoints();pts.SetData(numpy_to_vtk(xyz,deep=True));d=vtkPolyData();d.SetPoints(pts)
    pr=vtkProbeFilter();pr.SetInputData(d);pr.SetSourceData(data);pr.Update();d=pr.GetOutput()
    mask=vtk_to_numpy(d.GetPointData().GetArray('vtkValidPointMask')).astype(bool)
    u=vtk_to_numpy(d.GetPointData().GetArray('U')).copy();p=vtk_to_numpy(d.GetPointData().GetArray('p')).copy()
    assert np.isfinite(u).all() and np.isfinite(p).all()
    return mask,u,p

# Exactly the same deterministic equal-area surface points in both cases.
N=100000;i=np.arange(N);z=(i+.5)/N;theta=i*np.pi*(3-np.sqrt(5));rad=np.sqrt(1-z*z)
directions=np.column_stack([rad*np.cos(theta),rad*np.sin(theta),z])
all_summaries={}
for name in sys.argv[1:] or ['B00','A01']:
    C=R/name;D=O/name;D.mkdir(exist_ok=True);M=json.loads((C/'case-metadata.json').read_text());rho=M['rho_kg_m3']
    log=(C/'log.simpleFoam').read_text();assert '\nEnd\n' in log,'Solver has not finished normally'
    last=float(re.findall(r'^Time = (\S+)',log,re.M)[-1]);tn=f'{last:g}'
    check=(C/'log.checkMesh').read_text();extra=(C/'log.checkMesh-additional').read_text()
    assert 'Mesh OK' in check
    cells=int(re.search(r'cells:\s+(\d+)',check).group(1))
    h={}
    for lf in [C/'log.simpleFoam-warmup',C/'log.simpleFoam']:
        if not lf.exists():continue
        for block in re.split(r'\nTime = ',lf.read_text())[1:]:
            it=float(block.splitlines()[0]);vals={}
            for v,x in re.findall(r'Solving for (\w+), Initial residual = ([\deE+.-]+)',block):vals.setdefault(v,float(x))
            if all(v in vals for v in ['Ux','Uy','Uz','p','k','omega']):h[it]=dict(iteration=it,**{v:vals[v] for v in ['Ux','Uy','Uz','p','k','omega']})
    residuals=[h[k] for k in sorted(h)];csv_write(D/'residuals.csv',residuals)
    mon,port_areas=monitors(C);its=sorted(set.intersection(*[set(v) for v in mon.values()]))
    pp=M['ports'];hist=[]
    for it in its:
        ins=[p for p in pp if p['side']=='IN'];outs=[p for p in pp if p['side']=='OUT']
        hist.append(dict(iteration=it,inlet_mass_kg_s=-rho*sum(mon[p['patch']+'Flux'][it] for p in ins),outlet_mass_kg_s=rho*sum(mon[p['patch']+'Flux'][it] for p in outs),inlet_area_mean_static_Pa=rho*np.average([mon[p['patch']+'Pressure'][it] for p in ins],weights=[port_areas[p['patch']] for p in ins])))
    csv_write(D/'pressure-and-flow-history.csv',hist)
    boundary_history=[];regional_history=[]
    for it in its:
        for p in pp:
            patch=p['patch']
            boundary_history.append(dict(iteration=it,boundary=patch,region=p['region'],side=p['side'],native_area_m2=port_areas[patch],outward_mass_flow_kg_s=rho*mon[patch+'Flux'][it],area_mean_static_pressure_Pa=rho*mon[patch+'Pressure'][it]))
        for region in ['Water_Xminus','Water_Center','Water_Xplus']:
            ri=[p for p in pp if p['region']==region and p['side']=='IN'];ro=[p for p in pp if p['region']==region and p['side']=='OUT']
            qi=-rho*sum(mon[p['patch']+'Flux'][it] for p in ri);qo=np.array([rho*mon[p['patch']+'Flux'][it] for p in ro])
            pin=rho*np.average([mon[p['patch']+'Pressure'][it] for p in ri],weights=[port_areas[p['patch']] for p in ri])
            regional_history.append(dict(iteration=it,region=region,inlet_mass_kg_s=float(qi),outlet_mass_kg_s=float(qo.sum()),inlet_minus_outlet_kg_s=float(qi-qo.sum()),outlet_flow_spread_over_mean=float(np.ptp(qo)/qo.mean()),inlet_area_mean_static_Pa=float(pin)))
    csv_write(D/'boundary-monitor-history.csv',boundary_history)
    csv_write(D/'regional-pressure-and-flow-history.csv',regional_history)
    final_dir=C/'VTK'/f'{name}_{tn}'
    rr=vtkXMLUnstructuredGridReader();rr.SetFileName(str(final_dir/'internal.vtu'));rr.Update();data=rr.GetOutput()
    assert data.GetNumberOfCells()==cells
    U=vtk_to_numpy(data.GetCellData().GetArray('U'));assert np.isfinite(U).all()
    outrows=[];face_rows=[]
    for port in sorted([p for p in pp if p['side']=='OUT'],key=lambda p:p['patch']):
        patch=port['patch'];rd=vtkXMLPolyDataReader();rd.SetFileName(str(final_dir/'boundary'/f'{patch}.vtp'));rd.Update();d=rd.GetOutput()
        nf=d.GetNumberOfCells();assert nf>0
        ids=numpy_to_vtk(np.arange(nf),deep=True);ids.SetName('nativeFaceId');d.GetCellData().AddArray(ids)
        tri=vtkTriangleFilter();tri.SetInputData(d);tri.Update();t=tri.GetOutput()
        xyz=vtk_to_numpy(t.GetPoints().GetData());ff=vtk_to_numpy(t.GetPolys().GetConnectivityArray()).reshape(-1,3)
        av=np.cross(xyz[ff[:,1]]-xyz[ff[:,0]],xyz[ff[:,2]]-xyz[ff[:,0]])*.5
        assert (av[:,2]<0).all()
        face_ids=vtk_to_numpy(t.GetCellData().GetArray('nativeFaceId')).astype(int)
        area_vectors=np.column_stack([np.bincount(face_ids,weights=av[:,k],minlength=nf) for k in range(3)])
        a=np.linalg.norm(area_vectors,axis=1)
        assert (a>0).all()
        assert abs(a.sum()-port_areas[patch])/port_areas[patch]<2e-4,'VTK and native patch areas disagree'
        q=scalar_patch(C/tn/'phi',patch,nf)*rho;u=vtk_to_numpy(d.GetCellData().GetArray('U'));un=np.sum(u*area_vectors,axis=1)/a
        q_from_u=rho*np.sum(u*area_vectors,axis=1)
        boundary_p=vtk_to_numpy(d.GetCellData().GetArray('p'))
        assert np.max(np.abs(boundary_p))<1e-12,'Outlet gauge pressure not zero'
        if last in mon[patch+'Flux']:
            assert abs(q.sum()-rho*mon[patch+'Flux'][last])<1e-7,'Native flux and monitor disagree'
        for j in range(nf):face_rows.append(dict(outlet=patch,face_id=j,area_m2=a[j],native_mdot_kg_s=q[j],mdot_from_boundary_velocity_kg_s=q_from_u[j],outward_un_m_s=un[j]))
        outrows.append(dict(outlet=patch,region=port['region'],x_mm=port['x_mm'],y_mm=port['y_mm'],net_mass_kg_s=float(q.sum()),forward_mass_kg_s=float(q[q>0].sum()),reverse_mass_kg_s=float(-q[q<0].sum()),reverse_area_fraction=float(a[q<0].sum()/a.sum()),minimum_un_m_s=float(un.min()),mass_from_boundary_velocity_kg_s=float(q_from_u.sum()),flux_minus_boundary_velocity_mass_kg_s=float(q.sum()-q_from_u.sum())))
    csv_write(D/'18-outlets.csv',outrows);csv_write(D/'outlet-faces.csv',face_rows)
    older=sorted([p for p in (C/'processor0').iterdir() if re.fullmatch(r'\d+(\.\d+)?',p.name) and float(p.name)<last and (p/'U').exists()],key=lambda p:float(p.name))
    drift=[]
    for p in older[-2:]:
        num=den=0.
        for rank in range(4):
            a=native_u(C/f'processor{rank}'/p.name/'U');b=native_u(C/f'processor{rank}'/tn/'U');assert a.shape==b.shape
            num+=float(((a-b)**2).sum());den+=float((b*b).sum())
        drift.append(dict(from_iteration=float(p.name),to_iteration=last,unweighted_cell_vector_relative_L2=float(np.sqrt(num/den))))
    # Three interior surfaces. Near-wall surfaces remain distinct from no-slip walls.
    surfaces=[]
    prev_data=None;surface_iteration_stability=[]
    if older:
        prev_vtu=C/'VTK'/f'{name}_{older[-1].name}'/'internal.vtu'
        if prev_vtu.exists():
            prev_reader=vtkXMLUnstructuredGridReader();prev_reader.SetFileName(str(prev_vtu));prev_reader.Update();prev_data=prev_reader.GetOutput()
    for radius in [45.2,46.4,47.6]:
        xyz=directions*(radius*.001);valid,u,p=probe(data,xyz)
        speed=np.linalg.norm(u,axis=1);area_per_point=2*np.pi*radius**2/N
        stats=[]
        for threshold in [.1,.5,1.]:
            stats.append(dict(threshold_m_s=threshold,low_area_mm2=float((valid&(speed<threshold)).sum()*area_per_point),low_fraction_of_valid_fluid=float((speed[valid]<threshold).mean())))
        surfaces.append(dict(radius_mm=radius,valid_count=int(valid.sum()),invalid_count=int((~valid).sum()),area_mm2=float(valid.sum()*area_per_point),statistics=stats))
        np.savez_compressed(D/f'surface-r{radius:.1f}.npz',points_m=xyz,velocity_m_s=u,pressure_Pa=p*rho,valid=valid)
        if prev_data is not None:
            pm,pu,pprev=probe(prev_data,xyz);common=pm&valid
            sf=np.linalg.norm(u[common],axis=1);sp=np.linalg.norm(pu[common],axis=1)
            tests=[]
            for t in [.1,.5,1.]:
                af=sf<t;ap=sp<t;union=af|ap
                tests.append(dict(threshold_m_s=t,final_fraction=float(af.mean()),previous_fraction=float(ap.mean()),fraction_change_percentage_points=float(100*(af.mean()-ap.mean())),low_mask_intersection_over_union=float((af&ap).sum()/union.sum()) if union.any() else 1.,changed_mask_fraction=float((af!=ap).mean())))
            surface_iteration_stability.append(dict(radius_mm=radius,from_iteration=float(older[-1].name),to_iteration=last,common_points=int(common.sum()),relative_vector_L2_change=float(np.sqrt(((u[common]-pu[common])**2).sum()/(u[common]*u[common]).sum())),thresholds=tests))
    regions=[]
    for region in ['Water_Xminus','Water_Center','Water_Xplus']:
        q=np.array([r['net_mass_kg_s'] for r in outrows if r['region']==region]);assert len(q)==6
        ri=[p for p in pp if p['region']==region and p['side']=='IN']
        pin=rho*np.average([mon[p['patch']+'Pressure'][its[-1]] for p in ri],weights=[port_areas[p['patch']] for p in ri])
        regions.append(dict(region=region,total_outlet_kg_s=float(q.sum()),minimum_outlet_kg_s=float(q.min()),maximum_outlet_kg_s=float(q.max()),spread_over_mean=float(np.ptp(q)/q.mean()),CV=float(q.std()/q.mean()),inlet_area_mean_static_Pa=float(pin)))
    recent=[r for r in hist if r['iteration']>=last-200]
    summary=dict(name=name,actual_cells=cells,final_iteration=last,solver_ended_normally=True,SIMPLE_residual_goals_met='SIMPLE solution converged' in log,final_initial_residuals=residuals[-1],standard_mesh_check_passed=True,additional_mesh_check_passed='Mesh OK' in extra,final_inlet_mean_static_Pa=hist[-1]['inlet_area_mean_static_Pa'],recent_pressure_range_Pa=float(np.ptp([r['inlet_area_mean_static_Pa'] for r in recent])),final_net_mass_imbalance_kg_s=hist[-1]['inlet_mass_kg_s']-hist[-1]['outlet_mass_kg_s'],total_reverse_outlet_kg_s=sum(r['reverse_mass_kg_s'] for r in outrows),native_cell_velocity_iteration_changes=drift,midchannel_surfaces=surfaces,regional_outlet_distribution=regions,max_cell_speed_m_s=float(np.linalg.norm(U,axis=1).max()),CFD_model='Isothermal incompressible steady RANS k-omega SST',limitations=['No grid independence test or designed prism layers','Steady iteration variation is not physical time variation','Final snapshots cannot establish improvement if residuals/fields remain unsettled'])
    (D/'summary.json').write_text(json.dumps(summary,indent=2))
    summary['internal_surface_iteration_stability']=surface_iteration_stability
    summary['internal_iteration_monitor']=analyze_probe_history(C,D)
    summary['last_recorded_monitor_iteration']=hist[-1]['iteration']
    summary['native_port_area_m2']=port_areas
    summary['sum_abs_outlet_phi_vs_U_flux_difference_kg_s']=sum(abs(r['flux_minus_boundary_velocity_mass_kg_s']) for r in outrows)
    yp=re.findall(r'patch walls y\+\s*:\s*min\s*=\s*([\deE+.-]+),\s*max\s*=\s*([\deE+.-]+),\s*average\s*=\s*([\deE+.-]+)',log)
    summary['final_wall_yPlus']=dict(zip(['minimum','maximum','average'],map(float,yp[-1]))) if yp else None
    mv=re.search(r'Total volume\s*=\s*([\deE+.-]+)',check)
    summary['native_mesh_volume_m3']=float(mv.group(1).rstrip('.')) if mv else None
    cc=re.search(r'Concave cells[^\n]*number of cells:\s*(\d+)',extra)
    skew=re.search(r'Max skewness\s*=\s*([\deE+.-]+)',check)
    nonorth=re.search(r'non-orthogonality Max:\s*([\deE+.-]+)',check)
    summary['mesh_diagnostics']=dict(concave_cell_count=int(cc.group(1)) if cc else 0,maximum_skewness=float(skew.group(1)) if skew else None,maximum_nonorthogonality_degrees=float(nonorth.group(1)) if nonorth else None)
    (D/'summary.json').write_text(json.dumps(summary,indent=2))
    all_summaries[name]=summary
    print(json.dumps(summary,indent=2),flush=True)
    # Keep provenance and runnable solver settings, but no huge native volume mesh.
    import shutil
    for p in ['case-metadata.json','log.checkMesh','log.checkMesh-additional','log.simpleFoam','log.simpleFoam-warmup','log.potentialFoam']:
        if (C/p).exists():shutil.copy2(C/p,D/p)
    shutil.copytree(C/'system',D/'OpenFOAM-system',dirs_exist_ok=True)
    shutil.copytree(C/'constant',D/'OpenFOAM-properties',ignore=shutil.ignore_patterns('polyMesh','triSurface','extendedFeatureEdgeMesh'),dirs_exist_ok=True)
    initial=D/'OpenFOAM-initial-conditions';initial.mkdir(exist_ok=True)
    for field_name in ['U','p','k','omega','nut']:
        shutil.copy2(C/'0'/field_name,initial/field_name)
