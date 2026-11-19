"""Package only accepted, audited, visually reviewed scientific outputs."""
from pathlib import Path
import json,hashlib,zipfile,datetime
R=Path(__file__).parent
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
assessment=json.loads((R/'assessment.json').read_text());assert assessment['accepted']
assert json.loads((R/'figure_qa'/'visual_review.json').read_text())['passed']
exclude={'__pycache__','plot_cache','pdf_render','figure_drafts_v1','figure_drafts_v2','figure_drafts_v3','Pressure_bounds_pdf_preview.png','layout_preview.png','preview_layout.py','package_manifest.json','delivery_integrity.json'}
files=[p for p in sorted(R.rglob('*')) if p.is_file() and not any(part in exclude for part in p.relative_to(R).parts)]
manifest={'selected_run':assessment['selected_run'],'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'completed_runs':sorted(p.parent.name for p in (R/'runs').glob('*/summary.json')),
    'assessment':assessment,'files':[{'path':p.relative_to(R).as_posix(),'bytes':p.stat().st_size,'sha256':digest(p)} for p in files]}
(R/'package_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
files.append(R/'package_manifest.json')
deliveries={}
for name,chosen in [('Startup_corrected_simulations_and_figures.zip',files),
                    ('Figures_2-16_verified.zip',list(sorted((R/'figures').glob('*')))+[R/'Verification_report.txt',R/'Methods_for_manuscript.txt',R/'Final_convergence_table.csv'])]:
    target=R.parent/name
    if target.exists():raise FileExistsError('Refusing existing delivery '+str(target))
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p in chosen:z.write(p,p.relative_to(R).as_posix())
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        for p in chosen:
            h=hashlib.sha256()
            with z.open(p.relative_to(R).as_posix()) as f:
                for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
            assert h.hexdigest()==digest(p)
    deliveries[name]={'bytes':target.stat().st_size,'files':len(chosen),'sha256':digest(target),'CRC_and_all_file_hashes_verified':True}
    print(json.dumps({name:deliveries[name]}),flush=True)
(R/'delivery_integrity.json').write_text(json.dumps(deliveries,indent=2),encoding='utf-8')
(R.parent/'Startup_corrected_deliveries.sha256.txt').write_text('\n'.join(v['sha256']+'  '+k for k,v in deliveries.items())+'\n',encoding='utf-8')
