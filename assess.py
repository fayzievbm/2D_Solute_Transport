"""Require all completed comparisons and independent tests before figures."""
from pathlib import Path
import json,csv
R=Path(__file__).parent
def read(name):return json.loads((R/name).read_text()) if (R/name).exists() else {}
d=read('comparison_metrics.json')
required=['mesh0_dt2__domain128_dt2','mesh0_dt2__mesh1_dt2','mesh1_dt2__mesh2_dt2_t8',
          'mesh2_dt4_t8__mesh2_dt2_t8','mesh0_dt4__mesh0_dt2','mesh0_dt2__mesh0_dt1']
a={'accepted':False,'selected_run':'mesh2_dt2_t8','missing_comparisons':[k for k in required if k not in d],
   'criteria':{'domain_L2':.002,'time_L2':.002,'spatial_L2':.01,'spatial_integrals':.02,
   'spatial_point_over_peak':.02,'history_over_peak':.002,'mass_residual_fraction':1e-9,
   'pressure_water_balance_fraction':1e-6}}
if not a['missing_comparisons']:
    s=read('runs/mesh2_dt2_t8/summary.json');spatial=d[required[2]]
    histories=read('history_convergence.json');audits=read('output_audit.json')
    a['measured']={'domain_L2':d[required[0]]['maximum_field_l2'],
        'time_L2':max(d[k]['maximum_field_l2'] for k in required[3:]),
        'spatial_L2':spatial['maximum_field_l2'],
        'spatial_integrals':spatial['maximum_integral_change'],
        'spatial_point_over_peak':max(v['difference_over_field_peak'] for k,v in spatial['quantities'].items() if '_point_' in k),
        'history_over_peak':histories.get('maximum_difference_over_fine_peak',1.),
        'mass_residual_fraction':max(v['mass_residual_fraction'] for v in audits.values()) if audits else 1.,
        'pressure_water_balance_fraction':max(v['pressure_water_balance_fraction'] for v in audits.values()) if audits else 1.}
    a['criterion_passes']={k:a['measured'][k]<v for k,v in a['criteria'].items()}
    a['checks']={'decreasing_spatial_changes':d[required[2]]['maximum_field_l2']<d[required[1]]['maximum_field_l2'],
        'decreasing_temporal_changes_level0':d[required[5]]['maximum_field_l2']<d[required[4]]['maximum_field_l2'],
        'complete_histories':histories.get('complete',False),
        'independent_pressure_reference':read('verification/pressure_verification.json').get('passed',False),
        'independent_transport_reference':read('verification/transport_verification.json').get('passed',False),
        'pressure_tolerance_refinement':read('verification/pressure_tolerance_refinement.json').get('passed',False),
        'all_run_audits':all(audits.get(n,{}).get('checks_passed',False) for n in {n for pair in required for n in pair.split('__')}),
        'preservation':read('preservation_audit.json').get('sha256_verified',False)}
    a['accepted']=all(a['criterion_passes'].values()) and all(a['checks'].values())
    a['parameters']=s['parameters'];a['mesh']=s['mesh']
(R/'assessment.json').write_text(json.dumps(a,indent=2),encoding='utf-8')
with (R/'Convergence_summary.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['Coarse','Fine','Max_field_L2_percent','Inlet_region_L2_percent','Max_integral_change_percent'])
    for z in d.values():w.writerow([z['coarse'],z['fine'],100*z['maximum_field_l2'],100*z['maximum_inlet_region_l2'],100*z['maximum_integral_change']])
if a['accepted']:
    with (R/'Final_convergence_table.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['Check','Measured','Target_less_than','Unit','Result'])
        for k,v in a['measured'].items():
            scale=1 if k in ('mass_residual_fraction','pressure_water_balance_fraction') else 100
            w.writerow([k,scale*v,scale*a['criteria'][k],'fraction' if scale==1 else 'percent','PASS'])
print(json.dumps(a,indent=2))
