"""Independent matrix-exponential test of transport with known Darcy velocity."""
from pathlib import Path
from dataclasses import replace
import os,uuid,json
R=Path(__file__).parent;cache=R.parent/'scn'/uuid.uuid4().hex[:10];cache.mkdir(parents=True);os.environ['NUMBA_CACHE_DIR']=str(cache)
import numpy as np
from scipy.linalg import expm,solve
from corrected_solver import Model
from reference_solver import Parameters
from floating_point import flush_subnormals
from worker_floating_point import worker_subnormals
x=np.array([0,.01,.025,.05,.09,.14,.2,.28,.4]);y=np.linspace(0,.8,11)
p=replace(Parameters(),L=.4,inlet_width=.8,duration=2.,beta_ar=0.,beta_aa=0.,beta_ad=0.,beta_p0=0.,lambda_c=0.,lambda_a=0.,lambda_p=0.)
u=p.permeability_x*(p.p_inlet-p.p_ambient)/(p.viscosity*p.L*p.porosity)
nx=len(x)-2;ny=len(y)-2;hx=np.diff(x);hy=np.diff(y);vx=(hx[1:]+hx[:-1])/2;vy=(hy[1:]+hy[:-1])/2
xx,yy=np.meshgrid(x,y,indexing='ij');variable=p.D0*(1+np.exp(-p.omega_D*xx))*(1+4*yy/p.H*(1-yy/p.H))
exact={}
for case in ['constant','variable']:
    dx=(variable[:-1,1:-1]+variable[1:,1:-1])/2 if case=='variable' else np.full((nx+1,ny),p.D_x)
    dy=(variable[1:-1,:-1]+variable[1:-1,1:])/2 if case=='variable' else np.full((nx,ny+1),p.D_y)
    A=np.zeros((nx*ny,nx*ny));b=np.zeros(nx*ny)
    for i in range(nx):
        for j in range(ny):
            row=i*ny+j
            left=(u+dx[i,j]/hx[i])/vx[i];right=dx[i+1,j]/hx[i+1]/vx[i]
            A[row,row]-=left+right
            if i:A[row,row-ny]+=left
            else:b[row]+=left*p.c_inlet
            if i<nx-1:A[row,row+ny]+=right
            if j:
                lower=dy[i,j]/hy[j]/vy[j];A[row,row]-=lower;A[row,row-1]+=lower
            if j<ny-1:
                upper=dy[i,j+1]/hy[j+1]/vy[j];A[row,row]-=upper;A[row,row+1]+=upper
    ss=solve(A,-b);exact[case]=ss-expm(2*A)@ss
rows=[]
for dt in [.2,.1,.05]:
    m=Model(replace(p,dt=dt),x,y,all_cases=False)
    stationary=np.broadcast_to((p.p_inlet-p.p_ambient)*(1-x[1:-1,None]/p.L),(nx,ny)).copy()
    class FixedPressure:
        def at(self,t):return stationary.copy()
    m.pressure=FixedPressure();m.q=stationary.copy()
    with worker_subnormals(),flush_subnormals():
        for k in range(round(2/dt)):m.advance(dt)
    for index,case in enumerate(m.names):
        error=float(np.linalg.norm(m.states[index,0].ravel()-exact[case])/np.linalg.norm(exact[case]))
        mass=m.metrics(index)['relative_max_cumulative_residual'];assert mass<1e-10
        rows.append({'dt':dt,'case':case,'relative_L2_error':error,'mass_residual_fraction':mass})
for case in m.names:
    e=[r['relative_L2_error'] for r in rows if r['case']==case]
    assert all(a/b>3.5 for a,b in zip(e[:-1],e[1:])),(case,e)
result={'known_pore_velocity_m_per_s':u,'comparisons':rows,'passed':True,'interpretation':'Unchanged transport flux/reaction stages tested with zero reaction rates and a known steady Darcy field; independent full-matrix exponential provides the exact semidiscrete reference.'}
(R/'verification'/'transport_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
