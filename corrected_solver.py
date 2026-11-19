"""Implicit unsplit pressure plus conservative ADI transport on the fixed mesh."""
from pathlib import Path
from dataclasses import asdict,replace
import math,json,time,hashlib
import numpy as np
from blocked_solver import Model as State,coordinates
from reference_solver import Parameters
import blocked_transport as K
from pressure_implicit import Pressure

R=Path(__file__).parent

def startup_steps(dt,end=800.,initial=1e-4):
    growth=1+dt/end
    n=math.ceil(math.log1p(end*(growth-1)/initial)/math.log(growth))
    first=end*(growth-1)/math.expm1(n*math.log(growth))
    steps=first*np.exp(np.arange(n)*math.log(growth));steps*=end/float(steps.sum())
    assert dt/steps[-1]<=growth*(1+1e-6)
    return steps,growth

class Model(State):
    def __init__(self,p,x,y,all_cases=True,pressure_rtol=1e-7):
        super().__init__(p,x,y,all_cases)
        # Compile transport and flux signatures with a zero-source state.
        wx=np.zeros_like(self.pxa);wy=np.zeros_like(self.pya)
        K.transport(self.states[0,0],self.states[0,1],self.states[0,2],wx,wy,wy,
            self.Dx[0],self.Dy[0],self.hx,self.hy,self.vx,self.vy,self.f,1e-10,0.,self.r,self.stats[0])
        K.flux(self.q,self.pxa,self.pxb,0,p.p_inlet-p.p_ambient)
        K.flux(self.q,self.pya,self.pyb,1)
        self.pressure=Pressure(self,p.duration,rtol=pressure_rtol,atol=1e-7)
        self.time=0.;self.actual_steps=0

    def advance(self,dt):
        p=self.p;t=self.time;delta=dt/2
        qmid=self.pressure.at(t+delta);qn=self.pressure.at(t+dt)
        wx=p.storage*K.flux(qmid,self.pxa,self.pxb,0,p.p_inlet-p.p_ambient)
        wy0=p.storage*K.flux(self.q,self.pya,self.pyb,1)
        wy1=p.storage*K.flux(qn,self.pya,self.pyb,1)
        for s,phi in enumerate(self.phis):
            ih=p.c_inlet if phi==0 else p.c_inlet*(-np.expm1(-phi*(t+delta)))
            K.transport(self.states[s,0],self.states[s,1],self.states[s,2],wx,wy0,wy1,
                self.Dx[s],self.Dy[s],self.hx,self.hy,self.vx,self.vy,self.f,delta,ih,self.r,self.stats[s])
            it=p.c_inlet if phi==0 else p.c_inlet*(-np.expm1(-phi*t))
            for val in (it,ih):
                self.reservoir[s,0],self.reservoir[s,1],unused,unused2=K.react(val,self.reservoir[s,0],self.reservoir[s,1],delta,self.r)
        self.q=qn;self.time=t+dt;self.actual_steps+=1

def simulate(output,L=64.,level=2,dt=4.,duration=54000.,pressure_rtol=1e-7,all_cases=True):
    output=Path(output)
    if output.exists():raise FileExistsError('Refusing existing output: '+str(output))
    output.mkdir(parents=True)
    x,y=coordinates(L,level)
    p=replace(Parameters(),L=L,dx=.003125/2**level,dy=.0025/2**level,dt=dt,duration=duration,inlet_mode='fixed_width')
    start=time.perf_counter();m=Model(p,x,y,all_cases,pressure_rtol)
    print(output.name+': model constructed',flush=True)
    snapshots={name:{} for name in m.names};flow={};histories={name:[] for name in m.names}
    events=[t for t in (0.,1800.,3600.,7200.,18000.,54000.) if t<=duration]
    if duration not in events:events.append(duration)
    snapshot_steps={int(round(t/dt)):t for t in events}
    i=int(np.argmin(abs(x-.05)))-1;j=int(np.argmin(abs(y-.4)))-1
    assert abs(x[i+1]-.05)<1e-10 and abs(y[j+1]-.4)<1e-10
    next_event=0;last_history=-1e10;lastlog=time.perf_counter()
    def record(force=False):
        nonlocal next_event,last_history,lastlog
        t=m.time
        is_event=next_event<len(events) and abs(t-events[next_event])<1e-7
        if is_event:
            key=str(int(round(events[next_event]/dt)));flow[key]=dict(t=np.array(events[next_event]),**m.flow())
            for s,name in enumerate(m.names):snapshots[name][key]=dict(t=np.array(events[next_event]),**m.fields(s,events[next_event]))
            next_event+=1
        if force or is_event or t-last_history>=(1. if t<=800 else 8.)-1e-8:
            for s,name in enumerate(m.names):
                z=m.metrics(s);cc=m.states[s,:,i,j]
                histories[name].append([t,z['inventory'],z['cumulative_in'],z['cumulative_out'],z['cumulative_degradation'],z['final_residual'],z['min_c'],*cc])
            last_history=t
        now=time.perf_counter()
        if now-lastlog>=20:
            print(f'{output.name}: {t/3600:.4f}/{duration/3600:g} h; elapsed {now-start:.1f} s',flush=True);lastlog=now
    record(True)
    steps,growth=startup_steps(dt,initial=1e-4*(dt/4))
    from floating_point import flush_subnormals
    from worker_floating_point import worker_subnormals
    with worker_subnormals() as workers,flush_subnormals() as main:
        for h in steps:
            if m.time>=duration-1e-9:break
            m.advance(min(float(h),duration-m.time));record()
        if duration>=800:
            assert abs(m.time-800)<1e-7;m.time=800.
        while m.time<duration-1e-8:
            h=min(dt,duration-m.time)
            if next_event<len(events):h=min(h,events[next_event]-m.time)
            m.advance(h);record()
        record(True)
    common=dict(x=x,y=y,volume_x=m.vx,volume_y=m.vy,inlet_fraction=m.f)
    np.savez_compressed(output/'flow_snapshots.npz',**common,**{f'{k}_{field}':v for k,values in flow.items() for field,v in values.items()})
    for name in m.names:
        np.savez_compressed(output/(name+'_snapshots.npz'),**common,**{f'{k}_{field}':v for k,values in snapshots[name].items() for field,v in values.items()})
        # Deduplicate the final history row when it was already an event.
        history=histories[name]
        if len(history)>1 and history[-1][0]==history[-2][0]:history.pop()
        np.savetxt(output/(name+'_history.csv'),history,delimiter=',',header='time_s,inventory,cumulative_in,cumulative_out,cumulative_degradation,residual,min_c,c_at_005_mid,ca_at_005_mid,cp_at_005_mid',comments='')
    summary={'parameters':asdict(p),'scenario_names':m.names,'actual_transport_steps':m.actual_steps,'snapshot_steps':snapshot_steps,
        'snapshot_key_definition':'Requested physical time divided by subsequent transport dt; not the actual step count during startup.',
        'mesh':{'level':level,'nx':len(x),'ny':len(y),'min_dx':float(np.diff(x).min()),'max_dx':float(np.diff(x).max()),'min_dy':float(np.diff(y).min()),'max_dy':float(np.diff(y).max())},
        'effective_inlet_width':float(np.dot(m.f,m.vy)),'inventory_area':float(m.vx.sum()*m.vy.sum()),
        'time_integration':{'pressure':'Unsplit BDF/NDF, variable order 1-5, structured full-matrix linear solves','pressure_rtol':pressure_rtol,'pressure_atol_Pa':1e-7,
            'transport':'Conservative two-stage ADI with physical midpoint x velocity and endpoint y velocities; unchanged projected kinetic increments.',
            'startup_end_s':800.,'startup_growth':growth,'initial_step_s':float(steps[0]),'last_startup_step_s':float(steps[-1]),'subsequent_step_s':dt},
        'pressure':dict(m.pressure.stats,nlu=m.pressure.solver.nlu,nfev=m.pressure.solver.nfev),
        'transport':{name:m.metrics(s) for s,name in enumerate(m.names)},'floating_point':{'main':main,'workers':workers},
        'elapsed_seconds':time.perf_counter()-start,
        'code_hashes':{name:hashlib.sha256((R/name).read_bytes()).hexdigest() for name in ['corrected_solver.py','pressure_implicit.py','blocked_solver.py','blocked_transport.py','blocked_base.py','reference_solver.py','floating_point.py','worker_floating_point.py']}}
    (output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({'completed':output.name,'elapsed_s':summary['elapsed_seconds'],'pressure':summary['pressure'],'mass_residuals':{n:v['relative_max_cumulative_residual'] for n,v in summary['transport'].items()}}),flush=True)
    return summary
