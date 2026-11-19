"""Compare all four complete plotted histories, including startup."""
from pathlib import Path
import json
import numpy as np
R=Path(__file__).parent
out={'complete':False,'target_peak_normalized_difference':.002,'cases':{}}
coarse=R/'runs'/'mesh2_dt4_t8';fine=R/'runs'/'mesh2_dt2_t8'
if (coarse/'summary.json').exists() and (fine/'summary.json').exists():
    summary=json.loads((fine/'summary.json').read_text())
    for case in summary['scenario_names']:
        a=np.loadtxt(coarse/(case+'_history.csv'),delimiter=',',skiprows=1)
        b=np.loadtxt(fine/(case+'_history.csv'),delimiter=',',skiprows=1)
        rows={}
        for col,field in [(7,'c'),(8,'ca'),(9,'cp')]:
            ai=np.interp(b[:,0],a[:,0],a[:,col]);diff=abs(ai-b[:,col])
            peak=max(abs(b[:,col]).max(),1e-300);idx=int(np.argmax(diff))
            rows[field]={'maximum_difference_over_fine_peak':float(diff[idx]/peak),
                'maximum_absolute_difference':float(diff[idx]),'time_of_maximum_s':float(b[idx,0]),
                'rms_over_fine_peak':float(np.sqrt(np.mean(diff**2))/peak),
                'first_30min_difference_over_peak':float(diff[b[:,0]<=1800].max()/peak)}
        out['cases'][case]=rows
out['complete']=len(out['cases'])==4
if out['complete']:
    out['maximum_difference_over_fine_peak']=max(z['maximum_difference_over_fine_peak'] for rows in out['cases'].values() for z in rows.values())
    out['passed']=out['maximum_difference_over_fine_peak']<out['target_peak_normalized_difference']
out['interpretation']='Fixed physical point (0.05,0.4) m, 0-15 h. Coarse history linearly interpolated onto fine history. Startup samples are approximately 1 s apart; after 800 s both are sampled every 8 s, with additional exact figure-time samples. Includes sampling/interpolation effects. Does not establish all pointwise startup errors near the inlet.'
(R/'history_convergence.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
