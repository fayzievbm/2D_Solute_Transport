"""Short execution-speed/equivalence test; not a scientific production run."""
from pathlib import Path
import os,uuid,time,json
R=Path(__file__).parent
cache=R.parent/'scn'/uuid.uuid4().hex[:10];cache.mkdir(parents=True);os.environ['NUMBA_CACHE_DIR']=str(cache)
import numpy as np
from numba import set_num_threads
from blocked_solver import Model,coordinates
from reference_solver import Parameters
import blocked_transport as K
from analyze import interpolate
from floating_point import flush_subnormals
from worker_floating_point import worker_subnormals
x,y=coordinates(64,2);m=Model(Parameters(L=64),x,y,True)
with np.load(R/'verification'/'pressure_level2_tol1e-08.npz') as z:
    key=[k for k in z.files if k.startswith('54000')][0];q=z[key]
wx=m.p.storage*K.flux(q,m.pxa,m.pxb,0,400000.)
wy=m.p.storage*K.flux(q,m.pya,m.pyb,1)
for s,name in enumerate(m.names):
    with np.load(R/'runs'/'mesh0_dt2'/(name+'_snapshots.npz')) as z:
        for k,f in enumerate(['c','ca','cp']):m.states[s,k]=interpolate(z['27000_'+f],z['x'],z['y'],x[1:-1],y[1:-1])
initial=m.states.copy();rows=[];reference=None
for n in [1,2,4,8]:
    set_num_threads(n);m.states[:]=initial;m.stats[:]=0
    with worker_subnormals(),flush_subnormals():
        # Compile/warm independently of timing, then restore identical state.
        K.transport(*m.states[0],wx,wy,wy,m.Dx[0],m.Dy[0],m.hx,m.hy,m.vx,m.vy,m.f,1.,.01,m.r,m.stats[0])
        m.states[:]=initial;m.stats[:]=0;start=time.perf_counter()
        for k in range(12):
            for s in range(4):K.transport(*m.states[s],wx,wy,wy,m.Dx[s],m.Dy[s],m.hx,m.hy,m.vx,m.vy,m.f,1.,.01,m.r,m.stats[s])
        elapsed=time.perf_counter()-start
    if reference is None:reference=m.states.copy()
    row={'threads':n,'seconds_per_four_case_step':elapsed/12,'max_state_difference_to_single_thread':float(abs(m.states-reference).max())}
    rows.append(row);print(json.dumps(row),flush=True)
(R/'verification'/'thread_benchmark.json').write_text(json.dumps(rows,indent=2))
