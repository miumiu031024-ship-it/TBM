"""Read genuine OpenFOAM probes. Iteration history is NOT physical time."""
from pathlib import Path
import json,csv
import numpy as np

def analyze_probe_history(case:Path,output:Path):
    text=(case/'system/midchannelProbeLocations').read_text()
    xyz=np.fromstring(text.replace('(',' ').replace(')',' '),sep=' ').reshape(-1,3)
    n=len(xyz);frames={}
    for p in sorted((case/'postProcessing/midchannelMonitoringProbes').glob('*/U'),key=lambda p:float(p.parent.name)):
        for line in p.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith('#'):continue
            a=np.fromstring(line.replace('(',' ').replace(')',' '),sep=' ')
            if len(a)!=1+n*3:continue # A live writer may be midway through its last line.
            frames[float(a[0])]=a[1:].reshape(n,3)
    assert frames,'No complete probe frames'
    times=np.array(sorted(frames));u=np.stack([frames[t] for t in times])
    valid=np.isfinite(u).all(axis=2)&(abs(u)<1e100).all(axis=2)
    masks={'All':np.ones(n,dtype=bool),'Water_Xminus':xyz[:,0]<-.0232,'Water_Center':abs(xyz[:,0])<.0232,'Water_Xplus':xyz[:,0]>.0232}
    ap=2*np.pi*46.4**2/8192
    rows=[]
    for j,t in enumerate(times):
        for region,mask in masks.items():
            good=valid[j]&mask;s=np.linalg.norm(u[j,good],axis=1);assert len(s)>0
            row=dict(iteration=float(t),region=region,valid_points=int(good.sum()),mean_speed_m_s=float(s.mean()),vector_relative_L2_change_from_previous=None,previous_iteration=float(times[j-1]) if j else None)
            for threshold,key in [(.1,'0p1'),(.5,'0p5'),(1.,'1p0')]:
                row['low_fraction_'+key]=float((s<threshold).mean());row['low_area_mm2_'+key]=float((s<threshold).sum()*ap)
            if j:
                common=valid[j]&valid[j-1]&mask;a=u[j,common];b=u[j-1,common]
                den=float((a*a).sum());row['vector_relative_L2_change_from_previous']=float(np.sqrt(((a-b)**2).sum()/den)) if den else None
                af=np.linalg.norm(a,axis=1)<.5;bf=np.linalg.norm(b,axis=1)<.5;union=af|bf
                row['low_mask_0p5_intersection_over_union']=float((af&bf).sum()/union.sum()) if union.any() else 1.
            else:row['low_mask_0p5_intersection_over_union']=None
            rows.append(row)
    output.mkdir(parents=True,exist_ok=True)
    with (output/'internal-iteration-monitor.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    np.savez_compressed(output/'internal-iteration-monitor.npz',points_m=xyz,iterations=times,velocity_m_s=u,valid=valid)
    recent=[]
    for region in masks:
        rr=[q for q in rows if q['region']==region and q['iteration']>=times[-1]-100]
        changes=[q['vector_relative_L2_change_from_previous'] for q in rr if q['vector_relative_L2_change_from_previous'] is not None]
        recent.append(dict(region=region,frames=len(rr),from_iteration=rr[0]['iteration'],to_iteration=rr[-1]['iteration'],low_fraction_0p5_range_percentage_points=100*(max(q['low_fraction_0p5'] for q in rr)-min(q['low_fraction_0p5'] for q in rr)),maximum_successive_vector_relative_L2_change=max(changes) if changes else None))
    report=dict(original_equal_area_sample_count=8192,actual_clipped_point_count=n,radius_mm=46.4,divider_margin_mm=.01,first_recorded_iteration=float(times[0]),last_recorded_iteration=float(times[-1]),frames=len(times),recent_100_iteration_variation=recent,interpretation='Steady numerical iteration variation; NOT physical time fluctuations or time-averaged performance. Separate monitor sampling from the 100000-point final VTK comparison.')
    (output/'internal-iteration-monitor-summary.json').write_text(json.dumps(report,indent=2))
    return report

if __name__=='__main__':
    import sys
    root=Path(__file__).resolve().parent
    name=sys.argv[1]
    print(json.dumps(analyze_probe_history(root/name,root/'monitoring-live'/name),indent=2))
