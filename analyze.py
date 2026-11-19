"""Compare completed runs, with common physical-region weighted field norms."""
from pathlib import Path
import json,csv
import numpy as np
R=Path(__file__).parent;RUNS=R/'runs'
TIMES=[1800,3600,7200,18000,54000]

def load(run):
    s=json.loads((RUNS/run/'summary.json').read_text())
    z={name:np.load(RUNS/run/(name+'_snapshots.npz')) for name in s['scenario_names']}
    z['flow']=np.load(RUNS/run/'flow_snapshots.npz')
    return s,z

def interpolate(z,x,y,xx,yy):
    b=np.array([np.interp(xx,x,z[:,j]) for j in range(len(y))]).T
    return np.array([np.interp(yy,y,row) for row in b])

def compare(coarse,fine):
    a,az=load(coarse);b,bz=load(fine)
    x=bz['flow']['x'];y=bz['flow']['y'];ix=np.where((x>1e-12)&(x<=4+1e-12))[0];iy=np.arange(1,len(y)-1)
    xx=x[ix];yy=y[iy];vx=bz['flow']['volume_x'][ix-1].copy();vy=bz['flow']['volume_y']
    # The last comparison control volume is cut at the physical x=4 plane.
    if abs(xx[-1]-4)<1e-10:vx[-1]=(x[ix[-1]]-x[ix[-1]-1])/2
    w=vx[:,None]*vy[None,:];results={};small=xx<=.2+1e-12
    for case in [s for s in a['scenario_names'] if s in b['scenario_names']]+['flow']:
        for fld in (['p'] if case=='flow' else ['c','ca','cp']):
            for t in TIMES:
                ka=str(round(t/a['parameters']['dt']))+'_'+fld;kb=str(round(t/b['parameters']['dt']))+'_'+fld
                aa=az[case][ka];bb=bz[case][kb]
                # Both solutions evaluated at the finer solution's physical nodes.
                ai=interpolate(aa,az[case]['x'],az[case]['y'],xx,yy);bi=bb[np.ix_(ix,iy)]
                if fld=='p':ai=ai-a['parameters']['p_ambient'];bi=bi-b['parameters']['p_ambient']
                d=ai-bi;den=np.sum(w*bi*bi)
                z={'l2':float(np.sqrt(np.sum(w*d*d)/max(den,1e-300))),
                   'inlet_region_l2':float(np.sqrt(np.sum(w[small]*d[small]**2)/max(np.sum(w[small]*bi[small]**2),1e-300))),
                   'max_abs':float(abs(d).max()),'max_over_peak':float(abs(d).max()/max(abs(bi).max(),1e-300))}
                results[f'{case}_{fld}_{t}']=z
    quantities={}
    for case in [s for s in a['scenario_names'] if s in b['scenario_names']]:
        for f in ['inventory','cumulative_in','cumulative_degradation']:
            av=a['transport'][case][f];bv=b['transport'][case][f]
            quantities[f'{case}_{f}']={'coarse':av,'fine':bv,'relative_change':abs(av-bv)/max(abs(bv),1e-300)}
        for t in TIMES:
            for fld in ['c','ca','cp']:
                aa=az[case][f'{round(t/a["parameters"]["dt"])}_{fld}'];bb=bz[case][f'{round(t/b["parameters"]["dt"])}_{fld}']
                for xp,yp in [(0,0),(.05,.4),(.2,.4),(1,.4),(2,.4),(4,.4)]:
                    av=interpolate(aa,az[case]['x'],az[case]['y'],[xp],[yp])[0,0]
                    bv=interpolate(bb,bz[case]['x'],bz[case]['y'],[xp],[yp])[0,0]
                    quantities[f'{case}_{fld}_{t}_point_{xp}_{yp}']={'coarse':float(av),'fine':float(bv),'relative_change':float(abs(av-bv)/max(abs(bv),1e-300)),
                        'difference_over_field_peak':float(abs(av-bv)/max(abs(bb).max(),1e-300))}
    for z in [*az.values(),*bz.values()]:z.close()
    return {'coarse':coarse,'fine':fine,'fields':results,'quantities':quantities,
            'maximum_field_l2':max(v['l2'] for v in results.values()),
            'maximum_inlet_region_l2':max(v['inlet_region_l2'] for v in results.values()),
            'maximum_integral_change':max(v['relative_change'] for k,v in quantities.items() if '_point_' not in k)}

def main():
    # Only the new unsplit-pressure runs belong to this study.
    pairs=[('mesh0_dt2','domain128_dt2'),('mesh0_dt2','mesh1_dt2'),
           ('mesh1_dt2','mesh2_dt2_t8'),('mesh2_dt4_t8','mesh2_dt2_t8'),
           ('mesh0_dt4','mesh0_dt2'),('mesh0_dt2','mesh0_dt1')]
    results={}
    for a,b in pairs:
        if not (RUNS/a/'summary.json').exists() or not (RUNS/b/'summary.json').exists():continue
        z=compare(a,b);results[a+'__'+b]=z
        print(a,'->',b,'max L2 %',round(z['maximum_field_l2']*100,6),
              'inlet L2 %',round(z['maximum_inlet_region_l2']*100,6),'integral %',round(z['maximum_integral_change']*100,6),flush=True)
    (R/'comparison_metrics.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    with (R/'field_comparisons.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['Coarse','Fine','Field_time','L2_percent','Inlet_region_L2_percent','Max_over_peak_percent'])
        for name,z in results.items():
            for field,v in z['fields'].items():w.writerow([z['coarse'],z['fine'],field,100*v['l2'],100*v['inlet_region_l2'],100*v['max_over_peak']])
if __name__=='__main__':main()
