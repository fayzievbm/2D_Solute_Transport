"""Audit saved numerical fields, source provenance and preservation."""
from pathlib import Path
import json,hashlib
import numpy as np
R=Path(__file__).parent

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def audit():
    out={}
    for d in sorted((R/'runs').iterdir()):
        if not (d/'summary.json').exists():continue
        s=json.loads((d/'summary.json').read_text());p=s['parameters']
        assert s['scenario_names']==['constant','variable','exp_0.001','exp_0.0005']
        assert p['duration']==54000.
        with np.load(d/'flow_snapshots.npz') as z:
            x=z['x'];y=z['y'];vx=z['volume_x'];vy=z['volume_y'];f=z['inlet_fraction'];vol=vx[:,None]*vy[None,:]
            assert abs(float(np.dot(f,vy))-.04)<1e-12
            assert abs(x[-1]-p['L'])<1e-9 and abs(y[-1]-.8)<1e-12
            assert s['pressure']['min_gauge_pressure']>=-1e-6
            assert s['pressure']['max_gauge_pressure']<=p['p_inlet']-p['p_ambient']+1e-6
            max_velocity_error=0.;max_inventory_error=0.;max_cumulative=0.;max_saved_relative=0.
            for k,t in s['snapshot_steps'].items():
                pp=z[k+'_p'];wx=z[k+'_wx'];wy=z[k+'_wy']
                assert np.isfinite(pp).all() and np.isfinite(wx).all() and np.isfinite(wy).all()
                assert pp.min()>=p['p_ambient']-1e-6 and pp.max()<=p['p_inlet']+1e-6
                assert abs(pp[-1]-p['p_ambient']).max()<1e-7
                expectedx=-p['permeability_x']/p['viscosity']*np.diff(pp[:,1:-1],axis=0)/np.diff(x)[:,None]
                expectedy=-p['permeability_y']/p['viscosity']*np.diff(pp[1:-1,:],axis=1)/np.diff(y)[None,:]
                max_velocity_error=max(max_velocity_error,float(abs(wx-expectedx).max()),float(abs(wy-expectedy).max()))
                assert max_velocity_error<1e-12
                assert abs(pp[:,0]-pp[:,1]).max()<1e-8 and abs(pp[:,-1]-pp[:,-2]).max()<1e-8
        for case in s['scenario_names']:
            h=np.loadtxt(d/(case+'_history.csv'),delimiter=',',skiprows=1)
            assert np.isfinite(h).all() and np.all(np.diff(h[:,0])>0)
            nonzero=h[:,2]>0
            max_saved_relative=max(max_saved_relative,float(np.max(np.abs(h[nonzero,5])/h[nonzero,2])))
            with np.load(d/(case+'_snapshots.npz')) as z:
                for k,t in s['snapshot_steps'].items():
                    c,a,b=[z[k+'_'+v] for v in ['c','ca','cp']]
                    assert all(np.isfinite(v).all() for v in (c,a,b))
                    assert c.min()>=-1e-12 and a.min()>=0 and b.min()>=0
                    assert a.max()<=p['capacity_a'] and b.max()<=p['capacity_p']
                    inventory=float(np.sum(vol*(p['porosity']*c[1:-1,1:-1]+a[1:-1,1:-1]+b[1:-1,1:-1])))
                    row=h[np.argmin(abs(h[:,0]-t))];assert abs(row[0]-t)<1e-7
                    max_inventory_error=max(max_inventory_error,abs(inventory-row[1]))
            assert s['transport'][case]['min_c']>=-1e-12
            max_cumulative=max(max_cumulative,s['transport'][case]['relative_max_cumulative_residual'])
        assert max_inventory_error<1e-14 and max_cumulative<1e-9 and max_saved_relative<1e-9
        assert s['pressure']['max_relative_storage_integral_residual_after_1s']<1e-6
        for name,recorded_hash in s['code_hashes'].items():
            assert digest(R/name)==recorded_hash,(d.name,name,'source hash mismatch')
        out[d.name]={'mass_residual_fraction':max_cumulative,'snapshot_inventory_error':max_inventory_error,
            'darcy_velocity_error':max_velocity_error,'inlet_width':float(np.dot(f,vy)),
            'maximum_saved_residual_over_current_input':max_saved_relative,
            'pressure_water_balance_fraction':s['pressure']['max_relative_storage_integral_residual_after_1s'],
            'pressure_bounds_pass':True,'source_hashes_match':True,'checks_passed':True}
    (R/'output_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print('Completed runs audited:',len(out))

def preservation():
    previous=R.parent/'convergence_checks'
    manifest=json.loads((previous/'original_files_manifest.json').read_text())
    difference=json.loads((R/'preexisting_manifest_difference.json').read_text())
    errors=[name for name,v in manifest.items() if digest(name)!=v['sha256']]
    assert errors==[difference['path']],errors
    assert digest(difference['path'])==difference['current_sha256']
    results={'original_files_matching_previous_baseline':len(manifest)-1,
        'latest_Word_unchanged_since_current_digest_recorded':True,
        'historical_baseline_difference':difference,
        'sha256_verified':True}
    for folder in ['corrected_simulations','convergence_checks']:
        old=R.parent/folder;prior=json.loads((old/'package_manifest.json').read_text())
        errors=[v['path'] for v in prior['files'] if digest(old/v['path'])!=v['sha256']]
        assert not errors,(folder,errors)
        results[folder+'_files_unchanged']=len(prior['files'])
    (R/'preservation_audit.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results))

if __name__=='__main__':audit();preservation()
