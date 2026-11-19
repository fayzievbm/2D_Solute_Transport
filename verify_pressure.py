from pathlib import Path
from dataclasses import replace
import os,uuid,json,time
R=Path(__file__).parent
cache=R.parent/'scn'/uuid.uuid4().hex[:10];cache.mkdir(parents=True);os.environ['NUMBA_CACHE_DIR']=str(cache)
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import spsolve
from scipy.linalg import expm
from blocked_solver import Model
from reference_solver import Parameters
from pressure_implicit import Pressure,PressureOperator
from floating_point import flush_subnormals
import blocked_base as K
x=np.array([0,.005,.012,.025,.05,.1,.2,.4,1.])
y=np.array([0,.1,.2,.3,.35,.37,.38,.39,.4,.41,.42,.43,.45,.5,.6,.7,.8])
p=replace(Parameters(),L=1.,inlet_mode='fixed_width')
m=Model(p,x,y);o=PressureOperator(m);rng=np.random.default_rng(417)
q=rng.normal(size=m.q.shape)
direct=K.apply(q,m.pxa,m.pxb,0,m.vx,p.p_inlet-p.p_ambient)+K.apply(q,m.pya,m.pyb,1,m.vy)
err=float(np.linalg.norm(o.rhs(0,q.ravel())-direct.ravel())/np.linalg.norm(direct))
assert err<1e-13,err
solves={}
for gamma in [1e-6,.1,10.,10000.]:
    b=rng.normal(size=q.size);matrix=eye(q.size,format='csc')-gamma*o.A
    ref=spsolve(matrix,b);z=o.solve(o.factor(gamma),b)
    e=float(np.linalg.norm(z-ref)/np.linalg.norm(ref));assert e<1e-8,(gamma,e)
    solves[str(gamma)]={'relative_difference_to_sparse_LU':e,'relative_residual':float(np.linalg.norm(matrix@z-b)/np.linalg.norm(b))}
a=Pressure(m,100.,rtol=1e-9,atol=1e-9);b=Pressure(m,100.,rtol=1e-9,atol=1e-9,linear_solver='sparse')
steady=spsolve(o.A,-o.b);rows=[]
with flush_subnormals():
    for t in [1e-5,.001,.1,1.,10.,100.]:
        z=a.at(t).ravel();zz=b.at(t).ravel();exact=steady-expm(o.A.toarray()*t)@steady
        rel=float(np.linalg.norm(z-exact)/np.linalg.norm(exact));agreement=float(abs(z-zz).max())
        assert rel<2e-6,(t,rel);assert agreement<.01,(t,agreement)
        rows.append({'time_s':t,'relative_L2_to_matrix_exponential':rel,'max_Pa_difference_to_sparse_BDF':agreement})
assert a.stats['max_relative_storage_integral_residual_after_1s']<1e-6,a.stats
result={'operator_agreement':err,'linear_solves':solves,'time_integration':rows,'pressure_stats':a.stats,'passed':True}
(R/'verification'/'pressure_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
