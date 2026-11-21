"""Matched new OpenFOAM cases: baseline versus divider-root R0.25 mm.

STL and OpenFOAM units m; supplied CAD units mm. No previous CFD field reused.
"""
from pathlib import Path
import os, math, json, hashlib, importlib.util, shutil
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-a01-cache')
import cadquery as cq
import numpy as np
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepTools import BRepTools
from OCP.BRep import BRep_Tool
from OCP.TopLoc import TopLoc_Location
from OCP.TopAbs import TopAbs_REVERSED

ROOT=Path(__file__).resolve().parent
CAD=Path('/workspace/TBM/cad-proposals-20261120')
spec=importlib.util.spec_from_file_location('geometry',CAD/'generate_models.py')
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
rho=746.200955901462
nu=1.239443712658549e-7
base=g.read(CAD/'water/B00_baseline_water.step')

def facet(t):
    n=np.cross(t[1]-t[0],t[2]-t[0]);l=np.linalg.norm(n)
    if l<1e-25:return ''
    n/=l
    return 'facet normal '+' '.join(f'{x:.12g}' for x in n)+'\nouter loop\n'+''.join('vertex '+' '.join(f'{x:.12g}' for x in p)+'\n' for p in t)+'endloop\nendfacet\n'

def native_triangles(face):
    # Read the existing full-solid triangulation directly. CQ.tessellate can
    # re-mesh individual faces with RELATIVE deflection and reopen shared seams.
    loc=TopLoc_Location();poly=BRep_Tool.Triangulation_s(face.wrapped,loc)
    assert poly is not None
    v=np.array([poly.Node(i).Transformed(loc.Transformation()).Coord() for i in range(1,poly.NbNodes()+1)])*.001
    reverse=face.wrapped.Orientation()==TopAbs_REVERSED
    t=[]
    for i in range(1,poly.NbTriangles()+1):
        a,b,c=poly.Triangle(i).Get()
        t.append([a-1,c-1,b-1] if reverse else [a-1,b-1,c-1])
    return v,t

# An open ribbon around ORIGINAL divider/sphere intersection curves; used only
# for distance refinement, never as a wall or a snapping feature. Identical A/B.
guide=['solid rootGuide\n']
for solid in g.parts(base):
    for e in g.root_edges(solid):
        n=max(2,math.ceil(e.Length()/.10))
        v=np.array([e.positionAt(t).toTuple() for t in np.linspace(0,1,n)])*.001
        d=np.array([.000001,0,0])
        for a,b in zip(v[:-1],v[1:]):
            guide.extend([facet(np.array([a-d,a+d,b+d])),facet(np.array([a-d,b+d,b-d]))])
guide.append('endsolid rootGuide\n')
guide=''.join(guide)
monitor_n=8192
mi=np.arange(monitor_n);mz=(mi+.5)/monitor_n;mt=mi*np.pi*(3-np.sqrt(5));mr=np.sqrt(1-mz*mz)
monitor_xyz=np.column_stack((mr*np.cos(mt),mr*np.sin(mt),mz))*.0464
planes=[-23.699977677315708,-22.699977677315704,22.699999069163503,23.699999069163503]
mx=monitor_xyz[:,0]*1000
keep=np.ones(monitor_n,dtype=bool)
for lo,hi in zip(planes[::2],planes[1::2]):keep&=~((mx>lo-.01)&(mx<hi+.01))
monitor_xyz=monitor_xyz[keep]
monitor_locations='\n'.join('('+' '.join(f'{v:.14g}' for v in p)+')' for p in monitor_xyz)+'\n'
monitor_definition='''midchannelMonitoringProbes {type probes; libs (sampling); fields (U); interpolationScheme cellPoint;
executeControl timeStep; executeInterval 20; writeControl timeStep; writeInterval 20;
probeLocations (
#include "midchannelProbeLocations"
);}
'''

for name,file in [('B00','B00_baseline_water.step'),('A01','A01_root_R025_water.step')]:
    C=ROOT/name
    if C.exists():raise RuntimeError('Refuse overwrite '+str(C))
    (C/'constant/triSurface').mkdir(parents=True)
    (C/'system').mkdir();(C/'0').mkdir()
    (C/'system/midchannelProbeLocations').write_text(monitor_locations)
    s=g.read(CAD/'water'/file);ports=g.caps(s)
    assert len(ports)==36 and len(s.Solids())==3
    for p in ports:
        p['patch']=('inlet' if p['side']=='IN' else 'outlet')+p['port'][-2:]
    names=['walls',*[p['patch'] for p in ports]]
    streams={n:(C/'constant/triSurface'/f'{n}.stl').open('w') for n in [*names,'volume']}
    for n,f in streams.items():f.write('solid '+n+'\n')
    triangles=0
    for solid in g.parts(s):
        # Mesh the entire solid first so adjacent faces share the same edge
        # discretization. Meshing each face independently opens STL seams.
        BRepTools.Clean_s(solid.wrapped)
        BRepMesh_IncrementalMesh(solid.wrapped,.001,False,.12)
        for face in solid.Faces():
            patch='walls'
            if face.geomType()=='PLANE':
                for p in ports:
                    if np.linalg.norm(np.array(face.Center().toTuple())-[p[k] for k in ['x_mm','y_mm','z_mm']])<1e-5:
                        patch=p['patch'];break
            v,t=native_triangles(face)
            for tri in t:
                txt=facet(v[list(tri)]);streams[patch].write(txt);streams['volume'].write(txt);triangles+=1
    for n,f in streams.items():f.write('endsolid '+n+'\n');f.close()
    (C/'constant/triSurface/rootGuide.stl').write_text(guide)
    def write(p,text,cls='dictionary'):
        (C/p).write_text(f'FoamFile {{version 2.0; format ascii; class {cls}; object {Path(p).name};}}\n'+text+'\n')
    write('system/blockMeshDict','''scale 1;
vertices ((-.05 -.05 -.012) (.05 -.05 -.012) (.05 .05 -.012) (-.05 .05 -.012)
(-.05 -.05 .052) (.05 -.05 .052) (.05 .05 .052) (-.05 .05 .052));
blocks (hex (0 1 2 3 4 5 6 7) (28 28 18) simpleGrading (1 1 1));
edges (); boundary (background {type wall; faces ((0 3 2 1) (4 5 6 7) (0 1 5 4) (1 2 6 5) (2 3 7 6) (3 0 4 7));}); mergePatchPairs ();''')
    write('system/surfaceFeatureExtractDict','\n'.join(n+'.stl {extractionMethod extractFromSurface; extractFromSurfaceCoeffs {includedAngle 150;} writeObj no;}' for n in names))
    geometry='\n'.join(n+'.stl {type triSurfaceMesh; name '+n+';}' for n in [*names,'volume','rootGuide'])
    surfs='\n'.join(n+' {level (4 4); patchInfo {type '+('wall' if n=='walls' else 'patch')+';}}' for n in names)
    features=' '.join('{file "'+n+'.eMesh"; level 4;}' for n in names)
    seeds=[]
    for solid in g.parts(s):
        p=next(p for p in ports if p['region']==g.reg(solid) and p['side']=='IN')
        seed=[p['x_mm']*.001,p['y_mm']*.001,-.008];assert solid.isInside(tuple(np.array(seed)*1000));seeds.append(seed)
    locations=' '.join('(('+ ' '.join(map(str,p))+') '+n+')' for p,n in zip(seeds,['Water_Xminus','Water_Center','Water_Xplus']))
    write('system/snappyHexMeshDict',f'''castellatedMesh true; snap true; addLayers false;
geometry {{{geometry}}}
castellatedMeshControls {{maxLocalCells 3000000; maxGlobalCells 9000000; minRefinementCells 0; maxLoadUnbalance .1; nCellsBetweenLevels 3;
features ({features}); refinementSurfaces {{{surfs}}} resolveFeatureAngle 30;
refinementRegions {{volume {{mode inside; levels ((1e15 4));}} rootGuide {{mode distance; levels ((.00030 6) (.00060 5));}}}}
locationsInMesh ({locations}); allowFreeStandingZoneFaces true;}}
snapControls {{nSmoothPatch 5; tolerance 1.5; nSolveIter 60; nRelaxIter 8; nFeatureSnapIter 15; implicitFeatureSnap false; explicitFeatureSnap true; multiRegionFeatureSnap false;}}
addLayersControls {{relativeSizes true; layers {{}} expansionRatio 1.2; finalLayerThickness .3; minThickness .1; nGrow 0; featureAngle 60; nRelaxIter 5; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio .5; maxThicknessToMedialRatio .3; minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50;}}
meshQualityControls {{
#includeEtc "caseDicts/meshQualityDict"
minTetQuality 1e-9;
minVol 1e-18;
nSmoothScale 4; errorReduction .75;
}}
mergeTolerance 1e-7;''')
    write('system/decomposeParDict','numberOfSubdomains 4; method scotch;')
    write('constant/transportProperties',f'transportModel Newtonian; nu [0 2 -1 0 0 0 0] {nu:.16g};')
    write('constant/turbulenceProperties','simulationType RAS; RAS {RASModel kOmegaSST; turbulence on; printCoeffs on;}')
    def field(n,cls,dim,v,ib,ob,wb):
        write('0/'+n,f'dimensions {dim}; internalField uniform {v}; boundaryField {{\n'+'"inlet.*" {'+ib+'}\n"outlet.*" {'+ob+'}\nwalls {'+wb+'}\n}',cls)
    field('U','volVectorField','[0 1 -1 0 0 0 0]','(0 0 0)',f'type flowRateInletVelocity; massFlowRate constant {.23/18:.17g}; rhoInlet {rho:.17g}; value uniform (0 0 2.42);','type pressureInletOutletVelocity; value uniform (0 0 0);','type noSlip;')
    field('p','volScalarField','[0 2 -2 0 0 0 0]','0','type zeroGradient;','type fixedValue; value uniform 0;','type zeroGradient;')
    field('k','volScalarField','[0 2 -2 0 0 0 0]','.03','type turbulentIntensityKineticEnergyInlet; intensity .05; value uniform .03;','type inletOutlet; inletValue uniform .03; value uniform .03;','type kqRWallFunction; value uniform 1e-6;')
    field('omega','volScalarField','[0 0 -1 0 0 0 0]','1500','type turbulentMixingLengthFrequencyInlet; mixingLength .00021; value uniform 1500;','type inletOutlet; inletValue uniform 1500; value uniform 1500;','type omegaWallFunction; blending binomial; value uniform 1500;')
    field('nut','volScalarField','[0 2 -1 0 0 0 0]','0','type calculated; value uniform 0;','type calculated; value uniform 0;','type nutUSpaldingWallFunction; maxIter 100; tolerance 1e-7; value uniform 0;')
    write('system/fvSchemes','''ddtSchemes {default steadyState;} gradSchemes {default cellLimited Gauss linear 1;}
divSchemes {default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,k) bounded Gauss upwind; div(phi,omega) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear;}
laplacianSchemes {default Gauss linear limited .5;} interpolationSchemes {default linear;} snGradSchemes {default limited .5;}
wallDist {method meshWave;} fluxRequired {default no; p;}''')
    write('system/fvSolution','''solvers {p {solver GAMG; tolerance 1e-8; relTol .01; smoother DICGaussSeidel; cacheAgglomeration true; nCellsInCoarsestLevel 100; directSolveCoarsest yes;}
"(U|k|omega)" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol .05;}
Phi {solver GAMG; tolerance 1e-9; relTol 0; smoother GaussSeidel;}}
SIMPLE {nNonOrthogonalCorrectors 0; consistent yes; residualControl {p 1e-4; U 1e-5; k 1e-5; omega 1e-5;}}
potentialFlow {nNonOrthogonalCorrectors 5;}
relaxationFactors {fields {p .2;} equations {U .3; k .3; omega .3;}}''')
    funcs=[]
    for p in ports:
        for suffix,op,fld in [('Flux','sum','phi'),('Pressure','areaAverage','p')]:
            n=p['patch'];funcs.append(f'{n}{suffix} {{type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name {n}; operation {op}; fields ({fld}); writeFields false; log false; writeControl timeStep; writeInterval 10;}}')
    funcs.extend(['wallYPlus {type yPlus; libs (fieldFunctionObjects); executeControl writeTime; writeControl writeTime;}'])
    funcs.append(monitor_definition)
    write('system/controlDict','''application simpleFoam; startFrom startTime; startTime 0; stopAt endTime; endTime 500; deltaT 1;
writeControl timeStep; writeInterval 100; purgeWrite 3; writeFormat binary; writePrecision 8; writeCompression off; timeFormat general; timePrecision 8; runTimeModifiable true;
functions {'''+ '\n'.join(funcs)+'}')
    (C/'case.foam').touch()
    (C/'case-metadata.json').write_text(json.dumps(dict(name=name,source_STEP=file,source_sha256=hashlib.sha256((CAD/'water'/file).read_bytes()).hexdigest(),CAD_volume_mm3=g.volume(s),rho_kg_m3=rho,nu_m2_s=nu,mu_Pa_s=rho*nu,temperature_C=290,absolute_pressure_MPa=15.5,mass_flow_kg_s=.23,ports=ports,seeds_m=seeds,base_target_cells=5000000,STL_linear_tolerance_mm=.001,STL_relative_tolerance=False,nominal_root_cell_mm=.0558,nominal_bulk_cell_mm=.2232,explicit_prism_layers=False,triangles=triangles),indent=2))
    mf=C/'case-metadata.json';md=json.loads(mf.read_text())
    md['internal_iteration_monitor']=dict(radius_mm=46.4,original_equal_area_count=monitor_n,actual_fluid_point_count=len(monitor_xyz),divider_margin_mm=.01,write_every_iterations=20)
    md['model_assumptions']=dict(isothermal_constant_properties=True,gravity_included=False,smooth_walls=True)
    mf.write_text(json.dumps(md,indent=2))
    print(name,'prepared triangles',triangles,flush=True)
shutil.copy2(__file__,ROOT/'prepare-source-copy.py')
