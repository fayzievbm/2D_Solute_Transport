from pathlib import Path
from dataclasses import replace
import os,uuid,json,time,argparse
R=Path(__file__).parent
cache=R.parent/'scn'/uuid.uuid4().hex[:10];cache.mkdir(parents=True);os.environ['NUMBA_CACHE_DIR']=str(cache)
parser=argparse.ArgumentParser();parser.add_argument('--level',type=int,default=0);parser.add_argument('--rtol',type=float,default=1e-7);parser.add_argument('--until',type=float,default=54000.)
args=parser.parse_args()
import numpy as np
from blocked_solver import Model,coordinates
from reference_solver import Parameters
from pressure_implicit import Pressure
from floating_point import flush_subnormals
x,y=coordinates(64,args.level);m=Model(replace(Parameters(),L=64,inlet_mode='fixed_width'),x,y)
start=time.perf_counter();a=Pressure(m,args.until,rtol=args.rtol,atol=1e-7)
print('Constructed pressure system',m.q.shape,'seconds',time.perf_counter()-start,flush=True)
rows=[];snap={}
with flush_subnormals():
    for t in sorted(set([v for v in [1e-5,.001,.1,1,10,100,1000,1800,3600,7200,18000,54000] if v<=args.until]+[args.until])):
        q=a.at(t);rows.append({'time_s':t,'elapsed_s':time.perf_counter()-start,'steps':a.stats['adaptive_steps'],'nlu':a.solver.nlu,'max_pressure_Pa':float(q.max()+m.p.p_ambient)})
        if t>=1800:snap[str(int(t))]=q.copy()
        print(rows[-1],flush=True)
result={'level':args.level,'rtol':args.rtol,'rows':rows,'stats':a.stats,'passed':True}
name=f'pressure_level{args.level}_tol{args.rtol:g}'
(R/'verification'/(name+'.json')).write_text(json.dumps(result,indent=2),encoding='utf-8')
np.savez_compressed(R/'verification'/(name+'.npz'),x=x,y=y,**snap)
print(json.dumps(result['stats']),flush=True)
