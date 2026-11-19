"""Reproduce manuscript Table 4 from saved arrays; no simulation is rerun.

Usage: python normalized_inventory_analysis.py PATH_TO_STARTUP_CORRECTION NEW_OUTPUT_DIR
Requires numpy. All input files are opened read-only. Output directory must be new.
"""
from pathlib import Path
import csv,hashlib,json,sys
import numpy as np

def analyze(root):
    run=root/'runs/mesh2_dt2_t8'
    summary=json.loads((run/'summary.json').read_text())
    output={'concentration_convention':{
        'c':'C_water / C_ref','ca':'C_active_bulk / C_ref','cp':'C_passive_bulk / C_ref',
        'reference':'One common positive concentration scale; its physical value is unspecified.',
        'concentration_units':'dimensionless','coordinates':'m','time':'s',
        'integrated_M_I_O_G_units':'m^2',
        'dimensional_conversion':'C_ref times an integral gives parent-solute amount per unit transverse thickness.',
        'normalized_ratios':'dimensionless; percentages multiply by 100'},'cases':{},'input_hashes':{},'mesh_comparison':{}}
    for case in summary['scenario_names']:
        path=run/(case+'_snapshots.npz')
        output['input_hashes'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path) as z:
            volumes=z['volume_x'][:,None]*z['volume_y'][None,:]
            phase={name:float(np.sum(z['27000_'+name][1:-1,1:-1]*volumes)) for name in ('c','ca','cp')}
        phase['c']*=summary['parameters']['porosity']
        ledger=summary['transport'][case]
        with (run/(case+'_history.csv')).open(newline='') as fh:last=list(csv.DictReader(fh))[-1]
        assert float(last['time_s'])==54000
        for key in ('inventory','cumulative_in','cumulative_degradation'):
            assert abs(float(last[key])-ledger[key])<1e-16
        assert abs(sum(phase.values())-ledger['inventory'])<1e-14
        incoming=ledger['cumulative_in']
        percent={name:100*value/incoming for name,value in phase.items()}
        percent.update(stored=100*ledger['inventory']/incoming,degraded=100*ledger['cumulative_degradation']/incoming,outlet=100*ledger['cumulative_out']/incoming)
        assert abs(percent['stored']+percent['degraded']+percent['outlet']-100)<1e-8
        output['cases'][case]={'phase_integrals_m2':phase,'input_integral_m2':incoming,'percent_of_input':percent}
    for name in ('mesh0_dt2','mesh1_dt2','mesh2_dt2_t8','mesh2_dt4_t8'):
        path=root/'runs'/name/'summary.json'
        output['input_hashes'][name+'/summary.json']=hashlib.sha256(path.read_bytes()).hexdigest()
        data=json.loads(path.read_text())
        fractions={case:100*v['cumulative_degradation']/v['cumulative_in'] for case,v in data['transport'].items()}
        fractions['variable_minus_constant_percentage_points']=fractions['variable']-fractions['constant']
        output['mesh_comparison'][name]=fractions
    return output

if __name__=='__main__':
    if len(sys.argv)!=3:raise SystemExit(__doc__)
    source=Path(sys.argv[1]).resolve();dest=Path(sys.argv[2]).resolve()
    if dest.exists():raise SystemExit('Choose a new output directory; existing results are never overwritten.')
    data=analyze(source);dest.mkdir(parents=True)
    (dest/'normalized_inventory_results.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    with (dest/'Table_4_percentages.csv').open('w',newline='',encoding='utf-8') as fh:
        w=csv.writer(fh);w.writerow(['scenario','dissolved_percent','active_adsorbed_percent','passive_adsorbed_percent','stored_percent','degraded_percent'])
        for case,row in data['cases'].items():
            v=row['percent_of_input'];w.writerow([case]+[v[n] for n in ('c','ca','cp','stored','degraded')])
    print(json.dumps({'output_directory':str(dest),'checks_passed':True,'scenarios':list(data['cases'])},indent=2))
