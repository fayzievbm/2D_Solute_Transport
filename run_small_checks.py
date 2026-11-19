from pathlib import Path
import subprocess,sys
R=Path(__file__).parent
jobs=[('mesh0_dt2',0,64,2,1e-8),('domain128_dt2',0,128,2,1e-8),
      ('mesh0_dt4',0,64,4,1e-8),('mesh0_dt1',0,64,1,1e-8)]
for name,level,L,dt,tol in jobs:
    subprocess.run([sys.executable,str(R/'run.py'),'--name',name,'--level',str(level),'--L',str(L),
        '--dt',str(dt),'--pressure-rtol',str(tol)],cwd=R,check=True)
