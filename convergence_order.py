"""Three-grid observed orders/GCI for monotonic integral sequences only."""
from pathlib import Path
import json,math
R=Path(__file__).parent
results={}
families=[['mesh0_dt2','mesh1_dt2','mesh2_dt2_t8'],['mesh0_dt4','mesh0_dt2','mesh0_dt1']]
for names in families:
    if not all((R/'runs'/n/'summary.json').exists() for n in names):continue
    ss=[json.loads((R/'runs'/n/'summary.json').read_text()) for n in names];rows={}
    for case in [name for name in ss[0]['scenario_names'] if all(name in s['scenario_names'] for s in ss)]:
        for field in ['inventory','cumulative_in','cumulative_degradation']:
            vals=[s['transport'][case][field] for s in ss];d1=vals[0]-vals[1];d2=vals[1]-vals[2]
            z={'coarse_medium_fine':vals,'monotonic':d1*d2>0}
            if d1*d2>0 and abs(d1)>abs(d2)>0:
                p=math.log(abs(d1/d2))/math.log(2)
                z.update(observed_order=p,richardson_estimate=vals[2]+(vals[2]-vals[1])/(2**p-1),
                    fine_gci_fraction=1.25*abs((vals[2]-vals[1])/vals[2])/(2**p-1),
                    interpretation='Estimate conditional on asymptotic behaviour; not an independent proof of accuracy.')
            else:z['interpretation']='No positive monotonic observed order; no GCI claimed.'
            rows[case+'_'+field]=z
    results['__'.join(names)]=rows
(R/'observed_orders.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2))
