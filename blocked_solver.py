"""Same conservative ADI equations, fixed inlet, nested stretched far-field mesh.

This variant uses an inlet-focused mesh in both directions, with nested
subdivision of every interval. All spatial differences use actual face
distances and nodal control-volume widths. See coordinates() for the mesh.
"""
from pathlib import Path
from dataclasses import asdict,replace
import argparse,json,time,hashlib
import numpy as np
from reference_solver import Parameters
import blocked_transport as K

ROOT=Path(__file__).parent

def coordinates(L,level=0,tail_ratio=1.08,tail_cap=.8,uniform=False):
    if uniform:raise ValueError('This variant uses a nested inlet-focused mesh')
    def extend(values,end,h,ratio,cap):
        while values[-1]<end-1e-12:
            h=min(h*ratio,cap)
            values.append(min(values[-1]+h,end))
        return h
    x=list(np.linspace(0,.05,17))
    h=extend(x,.2,.003125,1.08,.0125)
    h=extend(x,2.,h,1.08,.025)
    h=extend(x,4.,h,1.08,.05)
    h=extend(x,L,h,tail_ratio,tail_cap)
    half=list(np.linspace(0,.05,21))
    extend(half,.4,.0025,1.12,.01)
    half=np.array(half);y=np.r_[.4-half[:0:-1],.4+half]
    def subdivide(a):
        n=2**level
        return np.r_[np.concatenate([np.linspace(u,v,n+1)[:-1] for u,v in zip(a[:-1],a[1:])]),a[-1]]
    return subdivide(np.array(x)),subdivide(y)

def full(z,f,inlet,outlet=0.):
    v=np.zeros((z.shape[0]+2,z.shape[1]+2));v[1:-1,1:-1]=z
    v[0,1:-1]=(1-f)*z[0,:]+f*inlet;v[-1,1:-1]=outlet
    v[:,0]=v[:,1];v[:,-1]=v[:,-2]
    return v

class Model:
    def __init__(self,p,x,y,all_cases=False):
        self.p=p;self.x=x;self.y=y
        self.hx=np.diff(x);self.hy=np.diff(y)
        self.vx=(self.hx[:-1]+self.hx[1:])/2
        self.vy=(self.hy[:-1]+self.hy[1:])/2
        lo=p.H/2-p.inlet_width/2;hi=p.H/2+p.inlet_width/2
        self.f=np.clip((np.minimum((y[1:-1]+y[2:])/2,hi)-np.maximum((y[:-2]+y[1:-1])/2,lo))/self.vy,0,1)
        nx=len(x)-2;ny=len(y)-2;xs=(nx+1,ny);ys=(nx,ny+1)
        self.q=np.zeros((nx,ny));self.names=['constant','variable']
        self.phis=np.array([0.,0.,.001,.0005] if all_cases else [0.,0.])
        if all_cases:self.names+=['exp_0.001','exp_0.0005']
        n=len(self.names);self.states=np.zeros((n,3,nx,ny));self.reservoir=np.zeros((n,2))
        self.stats=np.zeros((n,8));self.flowstats=np.zeros(3)
        self.pxa,self.pxb=K.coefficients(np.zeros(xs),np.full(xs,p.permeability_x/(p.viscosity*p.storage)),self.hx,0,self.f)
        self.pya,self.pyb=K.coefficients(np.zeros(ys),np.full(ys,p.permeability_y/(p.viscosity*p.storage)),self.hy,1,self.f)
        xx,yy=np.meshgrid(x,y,indexing='ij')
        D=p.D0*(1+np.exp(-p.omega_D*xx))*(1+4*yy/p.H*(1-yy/p.H))
        self.Dx=np.stack([(D[:-1,1:-1]+D[1:,1:-1])/2 if s=='variable' else np.full(xs,p.D_x) for s in self.names])
        self.Dy=np.stack([(D[1:-1,:-1]+D[1:-1,1:])/2 if s=='variable' else np.full(ys,p.D_y) for s in self.names])
        self.r=np.array([p.porosity,p.c_inlet,p.capacity_a,p.capacity_p,p.transition_a,p.transition_p,
                         p.beta_ar,p.beta_aa,p.beta_ad,p.beta_p0,p.lambda_c,p.lambda_a,p.lambda_p])
        self.step_index=0
    def step(self,n=1):
        p=self.p
        self.q=K.advance(self.q,self.states,self.reservoir,self.stats,self.flowstats,self.step_index,n,p.dt,
            self.hx,self.hy,self.vx,self.vy,self.f,self.pxa,self.pxb,self.pya,self.pyb,self.Dx,self.Dy,self.phis,
            self.r,p.p_inlet-p.p_ambient,p.storage)
        self.step_index+=n
    def metrics(self,s):
        p=self.p;c,a,b=self.states[s];z=self.stats[s]
        m=K.amount(c,a,b,self.vx,self.vy,p.porosity)
        return dict(inventory=m,cumulative_in=z[0],cumulative_out=z[1],cumulative_degradation=z[2],
            final_residual=m-z[0]+z[1]+z[2],max_cumulative_residual=z[3],
            relative_max_cumulative_residual=z[3]/max(z[0],1e-300),max_half_step_residual=z[4],
            min_c=z[5],max_c=z[6],capacity_projection_count=int(z[7]),max_ca=float(a.max()),max_cp=float(b.max()))
    def flow(self):
        p=self.p
        return dict(p=full(self.q,self.f,p.p_inlet-p.p_ambient)+p.p_ambient,
            wx=p.storage*K.flux(self.q,self.pxa,self.pxb,0,p.p_inlet-p.p_ambient),
            wy=p.storage*K.flux(self.q,self.pya,self.pyb,1))
    def fields(self,s,t):
        p=self.p;phi=self.phis[s];c0=p.c_inlet if phi==0 else p.c_inlet*(-np.expm1(-phi*t))
        c,a,b=self.states[s]
        return dict(c=full(c,self.f,c0),ca=full(a,self.f,self.reservoir[s,0]),cp=full(b,self.f,self.reservoir[s,1]))

def simulate(output,L=64.,level=0,dt=2.,duration=54000.,tail_ratio=1.08,tail_cap=.8,
             uniform=False,all_cases=False,snapshot_times=(0,1800,3600,7200,18000,54000)):
    output=Path(output)
    if output.exists():raise FileExistsError('Refusing to overwrite '+str(output))
    output.mkdir(parents=True)
    x,y=coordinates(L,level,tail_ratio,tail_cap,uniform)
    p=replace(Parameters(),L=L,dx=.003125/2**level,dy=.0025/2**level,dt=dt,duration=duration,inlet_mode='fixed_width')
    m=Model(p,x,y,all_cases);n=int(round(duration/dt));assert abs(n*dt-duration)<1e-8
    events={int(round(t/dt)):t for t in snapshot_times if t<=duration};events[n]=duration
    snapshots={s:{} for s in m.names};flow={};history={s:[] for s in m.names}
    point_i=int(np.argmin(abs(x-.05)))-1;point_j=int(np.argmin(abs(y-.4)))-1
    assert abs(x[point_i+1]-.05)<1e-10
    start=time.perf_counter();lastlog=start
    while True:
        k=m.step_index;t=k*dt
        if k in events:
            flow[str(k)]=dict(t=np.array(t),**m.flow())
            for s,name in enumerate(m.names):snapshots[name][str(k)]=dict(t=np.array(t),**m.fields(s,t))
        for s,name in enumerate(m.names):
            mm=m.metrics(s);cc=m.states[s,:,point_i,point_j]
            history[name].append([t,mm['inventory'],mm['cumulative_in'],mm['cumulative_out'],mm['cumulative_degradation'],mm['final_residual'],mm['min_c'],*cc])
        if k>=n:break
        todo=min(max(1,int(round(10/dt))),n-k)
        future=[a-k for a in events if a>k]
        if future:todo=min(todo,min(future))
        m.step(todo)
        now=time.perf_counter()
        if now-lastlog>20 or m.step_index==n:
            print(f'{output.name}: {m.step_index*dt/3600:.3f}/15 h, elapsed {now-start:.1f} s',flush=True);lastlog=now
    common=dict(x=x,y=y,volume_x=m.vx,volume_y=m.vy,inlet_fraction=m.f)
    np.savez_compressed(output/'flow_snapshots.npz',**common,**{f'{k}_{f}':v for k,dd in flow.items() for f,v in dd.items()})
    for name in m.names:
        np.savez_compressed(output/(name+'_snapshots.npz'),**common,**{f'{k}_{f}':v for k,dd in snapshots[name].items() for f,v in dd.items()})
        np.savetxt(output/(name+'_history.csv'),history[name],delimiter=',',header='time_s,inventory,cumulative_in,cumulative_out,cumulative_degradation,residual,min_c,c_at_005_mid,ca_at_005_mid,cp_at_005_mid',comments='')
    summary=dict(parameters=asdict(p),scenario_names=m.names,steps=n,snapshot_steps=events,
        mesh=dict(inlet_focused=True,level=level,max_dy=float(np.diff(y).max()),tail_ratio=tail_ratio,tail_cap=tail_cap,uniform=uniform,nx=len(x),ny=len(y),
                  near_dx=p.dx,dy=p.dy,max_dx=float(np.diff(x).max())),
        effective_inlet_width=float(np.dot(m.f,m.vy)),inventory_area=float(m.vx.sum()*m.vy.sum()),
        inventory_definition='sum(volume_x*volume_y*(theta*c+ca+cp)) over interior nodal control volumes',
        flow=dict(max_storage_flux_residual=float(m.flowstats[0]),min_pressure=p.p_ambient+float(m.flowstats[1]),max_pressure=p.p_ambient+float(m.flowstats[2])),
        transport={name:m.metrics(s) for s,name in enumerate(m.names)},elapsed_seconds=time.perf_counter()-start,
        code_hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['blocked_base.py','blocked_transport.py','blocked_solver.py','reference_solver.py']})
    (output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({'run':output.name,'elapsed':summary['elapsed_seconds'],'mass_residuals':{k:v['relative_max_cumulative_residual'] for k,v in summary['transport'].items()}}),flush=True)
    return summary

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--name',required=True);a.add_argument('--L',type=float,default=64)
    a.add_argument('--level',type=int,default=0);a.add_argument('--dt',type=float,default=2);a.add_argument('--duration',type=float,default=54000)
    a.add_argument('--tail-ratio',type=float,default=1.08);a.add_argument('--tail-cap',type=float,default=.8)
    a.add_argument('--uniform',action='store_true');a.add_argument('--all-cases',action='store_true')
    args=vars(a.parse_args());name=args.pop('name');simulate(ROOT/'runs'/name,**args)
