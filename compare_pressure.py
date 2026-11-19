from pathlib import Path
import json,numpy as np
R=Path(__file__).parent;rows=[]
with np.load(R/'verification'/'pressure_level2_tol1e-07.npz') as a,np.load(R/'verification'/'pressure_level2_tol1e-08.npz') as b:
    x=b['x'];y=b['y'];vx=(np.diff(x)[:-1]+np.diff(x)[1:])/2;vy=(np.diff(y)[:-1]+np.diff(y)[1:])/2
    select=x[1:-1]<=4+1e-10;weights=vx[select,None]*vy[None,:]
    for t in ['1800','3600','7200','18000','54000']:
        aa=a[t][select];bb=b[t][select];err=aa-bb
        rows.append({'time_s':int(t),'weighted_relative_L2':float(np.sqrt(np.sum(weights*err**2)/np.sum(weights*bb**2))),'max_absolute_difference_Pa':float(abs(err).max())})
out={'comparisons':rows,'max_field_L2':max(r['weighted_relative_L2'] for r in rows),'target':2e-3}
out['passed']=out['max_field_L2']<out['target']
(R/'verification'/'pressure_tolerance_refinement.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
