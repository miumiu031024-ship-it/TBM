"""Long first-order startup then second-order restoration; no geometry change."""
from pathlib import Path
import shutil,json,re,hashlib
R=Path(__file__).resolve().parent;S=Path('/workspace/scratch/tbm-cause-20261121/S30_upwind');C=R/'F01_upwind_long';assert not C.exists()
(C/'constant').mkdir(parents=True);(C/'600').mkdir();(C/'system').mkdir()
mesh=Path('/workspace/scratch/tbm-cause-20261121/SeedSide/constant/polyMesh');(C/'constant/polyMesh').symlink_to(mesh,target_is_directory=True)
for n in ['transportProperties','turbulenceProperties']:shutil.copy2(S/'constant'/n,C/'constant'/n)
for n in ['U','p','phi','k','omega','nut']:shutil.copy2(S/'600'/n,C/'600'/n)
for n in ['fvSchemes','fvSolution','decomposeParDict']:shutil.copy2(S/'system'/n,C/'system'/n)
s=(S/'system/controlDict').read_text().replace('startTime 500;','startTime 600;').replace('endTime 600;','endTime 800;').replace('writeInterval 50; purgeWrite 0;','writeInterval 100; purgeWrite 3;')
# Pure monitoring changes; no field updates or physical-model changes.
s=s.replace('executeInterval 5;','executeInterval 20;').replace('writeInterval 5;','writeInterval 20;').replace('writeInterval 1;}','writeInterval 5;}')
(C/'system/controlDict').write_text(s)
meta=dict(parent=str(S),selected_region='Water_Xminus',cells=1264925,original_total_inlet_kg_s=.23,selected_inlet_kg_s=.23/3,inlets=6,outlets=6,temperature_C=290,absolute_reference_MPa=15.5,rho_kg_m3=746.200955901462,nu_m2_s=1.239443712658549e-7,first_order_start_iteration=600,geometry_and_mesh_changed=False,steady_iterations_are_not_physical_time=True)
meta['ports']=[p for p in json.loads(Path('/workspace/scratch/tbm-cause-20261121/SeedSide/case-metadata.json').read_text())['ports'] if p['region']=='Water_Xminus']
meta['restart_field_sha256']={n:hashlib.sha256((C/'600'/n).read_bytes()).hexdigest() for n in ['U','p','phi','k','omega','nut']}
meta['mesh_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in mesh.iterdir() if p.is_file()}
(C/'case-metadata.json').write_text(json.dumps(meta,indent=2))
plan=dict(first_order=dict(case='F01_upwind_long',start_iteration=600,checkpoint_interval=200,minimum_extra_iterations=400,maximum_iteration=2000,stop_on_stall=False),second_order=dict(case='F02_second_order_restore',maximum_extra_iterations=600,checkpoint_interval=200),assessment=dict(residual_max_last100=dict(p=1e-4,Ux=1e-5,Uy=1e-5,Uz=1e-5,k=1e-5,omega=1e-5),native_velocity_relative_L2_last100_below=.01,maximum_20_iteration_probe_vector_relative_L2_last100_below=.01,probe_low_fraction_range_last100_pp_below=.2,each_outlet_flow_range_last100_fraction_below=.001,absolute_total_mass_imbalance_fraction_below=.001,all_required=True),claims_excluded=['Physical time/frequency from steady iterations','Grid independence','Thermal or manufacturing validation','STAR raw-data reproduction','A01 performance without matched verification'])
(R/'plan.json').write_text(json.dumps(plan,indent=2))
print('Prepared unchanged side-domain first-order continuation at iteration600')
