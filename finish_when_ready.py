"""Current-task pipeline; waits for planned runs, then applies acceptance gates."""
from pathlib import Path
import time,json,subprocess,sys
R=Path(__file__).parent
required=['mesh0_dt2','domain128_dt2','mesh0_dt4','mesh0_dt1','mesh1_dt2','mesh2_dt4_t8','mesh2_dt2_t8']
last=None
while True:
    ready=[]
    for name in required:
        try:
            s=json.loads((R/'runs'/name/'summary.json').read_text())
            if s['parameters']['duration']==54000:ready.append(name)
        except (FileNotFoundError,json.JSONDecodeError):pass
    if ready!=last:print('Completed planned runs: '+str(len(ready))+'/7; remaining '+', '.join(n for n in required if n not in ready),flush=True);last=ready
    if len(ready)==len(required):break
    time.sleep(10)
for script in ['analyze.py','convergence_order.py','assess_histories.py','audit_results.py','assess.py']:
    print('Running '+script,flush=True)
    subprocess.run([sys.executable,str(R/script)],cwd=R,check=True)
assert json.loads((R/'assessment.json').read_text())['accepted'],'Numerical acceptance failed; no final figures generated.'
for script in ['make_verified_figures.py','qa_figures.py','make_report.py']:
    print('Running '+script,flush=True)
    subprocess.run([sys.executable,str(R/script)],cwd=R,check=True)
print('Numerical gates and export checks passed. Manual visual review is still required before packaging.',flush=True)
