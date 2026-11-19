"""Conservative implementation of revised manuscript equations (1)-(53).

Original files are never imported or written. Units: m, s, Pa; concentrations
share the source model's amount basis. The state is on interior grid nodes.
Boundary nodes are algebraic boundary values; the discrete inventory is over
the solved control volumes, each of area dx*dy. See README.txt for scope.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import json, math, time
import numpy as np
from scipy.linalg.lapack import dgtsv, dgttrf, dgttrs

@dataclass(frozen=True)
class Parameters:
    L: float = 4.
    H: float = .8
    dx: float = .05
    dy: float = .04
    dt: float = 2.                    # full step; each directional stage dt/2
    duration: float = 54000.
    p_ambient: float = 1e5
    p_inlet: float = 5e5
    permeability_x: float = 5e-14
    permeability_y: float = 1e-14
    viscosity: float = 1e-3
    storage: float = 1e-8
    porosity: float = .3
    c_inlet: float = .01
    D_x: float = 5e-5
    D_y: float = 1e-5
    D0: float = 1e-5
    omega_D: float = 1.
    capacity_a: float = .004
    capacity_p: float = .008
    transition_a: float = .0002
    transition_p: float = .003
    beta_ar: float = 3e-5
    beta_aa: float = 2.5e-5
    beta_ad: float = 3e-5
    beta_p0: float = 7.5e-5            # both passive branches; source betp0
    lambda_c: float = 1e-5
    lambda_a: float = 2e-5
    lambda_p: float = 1.5e-5
    inlet_mode: str = "single_node"  # manuscript main run
    inlet_width: float = .04        # only for explicitly labelled sensitivity

class Grid:
    def __init__(self, p):
        self.p=p
        self.I=int(round(p.L/p.dx)); self.J=int(round(p.H/p.dy))
        assert abs(self.I*p.dx-p.L)<1e-10 and abs(self.J*p.dy-p.H)<1e-10
        assert self.J%2==0
        self.shape=(self.I-1,self.J-1)
        self.x=np.arange(self.I+1)*p.dx;self.y=np.arange(self.J+1)*p.dy
        self.cell_area=p.dx*p.dy
        self.inlet_fraction=np.zeros(self.J-1)
        if p.inlet_mode=="single_node":
            self.inlet_fraction[self.J//2-1]=1.
        elif p.inlet_mode=="fixed_width":
            y=self.y[1:-1]; lo=p.H/2-p.inlet_width/2;hi=p.H/2+p.inlet_width/2
            self.inlet_fraction=np.clip((np.minimum(y+p.dy/2,hi)-np.maximum(y-p.dy/2,lo))/p.dy,0,1)
        elif p.inlet_mode=="all": self.inlet_fraction[:]=1.
        else: raise ValueError(p.inlet_mode)
    def full(self, interior, inlet, outlet=0.):
        z=np.zeros((self.I+1,self.J+1))
        z[1:-1,1:-1]=interior
        f=self.inlet_fraction
        z[0,1:-1]=(1-f)*interior[0,:]+f*inlet
        z[-1,1:-1]=outlet
        z[:,0]=z[:,1];z[:,-1]=z[:,-2]
        return z

class Direction:
    """Lz = -div(F), F = a*z_left+b*z_right; single shared face flux."""
    def __init__(self, grid, axis, velocity, D):
        self.grid=grid;self.axis=axis
        h=grid.p.dx if axis==0 else grid.p.dy
        self.h=h
        self.a=np.maximum(velocity,0)+D/h
        self.b=np.minimum(velocity,0)-D/h
        if axis==0:
            f=grid.inlet_fraction
            # Velocity already includes open-face fraction; only dispersive
            # conductance is multiplied by fraction at a partial inlet face.
            self.a[0,:]=np.maximum(velocity[0,:],0)+D[0,:]/h*f
            self.b[0,:]=np.minimum(velocity[0,:],0)-D[0,:]/h*f
        else:
            self.a[:,0]=self.a[:,-1]=0
            self.b[:,0]=self.b[:,-1]=0
        if axis==0:
            self.lower=self.a[:-1,:]/h
            self.upper=-self.b[1:,:]/h
            self.diagonal=(self.b[:-1,:]-self.a[1:,:])/h
        else:
            self.lower=self.a[:,:-1]/h
            self.upper=-self.b[:,1:]/h
            self.diagonal=(self.b[:,:-1]-self.a[:,1:])/h
    def flux(self,z,left=0.,right=0.):
        if self.axis==0:
            f=np.empty_like(self.a)
            f[0,:]=self.a[0,:]*left+self.b[0,:]*z[0,:]
            f[-1,:]=self.a[-1,:]*z[-1,:]+self.b[-1,:]*right
            f[1:-1,:]=self.a[1:-1,:]*z[:-1,:]+self.b[1:-1,:]*z[1:,:]
        else:
            f=np.zeros_like(self.a)
            f[:,1:-1]=self.a[:,1:-1]*z[:,:-1]+self.b[:,1:-1]*z[:,1:]
        return f
    def apply(self,z,left=0.,right=0.):
        return -np.diff(self.flux(z,left,right),axis=self.axis)/self.h
    def source(self,left=0.,right=0.):
        s=np.zeros(self.grid.shape)
        if self.axis==0:
            s[0,:]+=self.a[0,:]*left/self.h
            s[-1,:]-=self.b[-1,:]*right/self.h
        return s
    def pack(self,x):
        return np.asarray(x).ravel(order="F" if self.axis==0 else "C").copy()
    def unpack(self,x):
        return np.asarray(x).reshape(self.grid.shape,order="F" if self.axis==0 else "C")
    def matrix(self,delta,decay=0.):
        lo=-delta*self.lower.copy();up=-delta*self.upper.copy()
        if self.axis==0:lo[0,:]=0;up[-1,:]=0
        else:lo[:,0]=0;up[:,-1]=0
        return self.pack(lo)[1:],self.pack(1+delta*decay-delta*self.diagonal),self.pack(up)[:-1]
    def factor(self,delta,decay=0.):
        lo,d,up=self.matrix(delta,decay)
        self.lu=dgttrf(lo,d,up)
        if self.lu[-1]!=0:raise ArithmeticError("Pressure factorization failed")
    def solve(self,rhs,delta,decay=0.,factored=False):
        b=self.pack(rhs).reshape(-1,1)
        if factored:
            x,info=dgttrs(*self.lu[:-1],b,overwrite_b=1)
        else:
            lo,d,up=self.matrix(delta,decay)
            _,_,_,x,info=dgtsv(lo,d,up,b,overwrite_dl=1,overwrite_d=1,overwrite_du=1,overwrite_b=1)
        if info!=0:raise ArithmeticError(f"Tridiagonal solve failed: {info}")
        return self.unpack(x[:,0])

class Pressure:
    def __init__(self,g):
        self.g=g;p=g.p
        self.delta=p.dt/2
        self.q=np.zeros(g.shape)   # gauge pressure, avoids cancellation of p0
        self.x=Direction(g,0,np.zeros((g.I,g.J-1)),np.full((g.I,g.J-1),p.permeability_x/(p.viscosity*p.storage)))
        self.y=Direction(g,1,np.zeros((g.I-1,g.J)),np.full((g.I-1,g.J),p.permeability_y/(p.viscosity*p.storage)))
        self.x.factor(self.delta);self.y.factor(self.delta)
        self.jump=p.p_inlet-p.p_ambient
    def velocities(self,q):
        # Pressure-operator flux is -(k/(mu*Sp))*grad p. Multiplication by Sp
        # yields Darcy velocity, including each wall/inlet face convention.
        return self.g.p.storage*self.x.flux(q,self.jump),self.g.p.storage*self.y.flux(q)
    def step(self):
        d=self.delta;q0=self.q
        rhs=q0+d*self.y.apply(q0)+d*self.x.source(self.jump)
        qs=self.x.solve(rhs,d,factored=True)
        q1=self.y.solve(qs+d*self.x.apply(qs,self.jump),d,factored=True)
        wx,unused=self.velocities(qs)
        unused,wy0=self.velocities(q0)
        unused,wy1=self.velocities(q1)
        # Independent storage/physical boundary-flux ledger for pressure.
        s=self.g.p.storage*self.g.cell_area
        storage_change=s*float(np.sum(q1-q0))
        flow_in=pore_boundary(wx,wy0,self.g)[0]+pore_boundary(wx,wy1,self.g)[0]
        flow_out=pore_boundary(wx,wy0,self.g)[1]+pore_boundary(wx,wy1,self.g)[1]
        residual=storage_change-d*(flow_in-flow_out)
        self.q=q1
        return wx,wy0,wy1,residual

def pore_boundary(fx,fy,g):
    oriented=np.concatenate((fx[0,:]*g.p.dy,-fx[-1,:]*g.p.dy,
                             fy[:,0]*g.p.dx,-fy[:,-1]*g.p.dx))
    return float(np.maximum(oriented,0).sum()),float(-np.minimum(oriented,0).sum())

def reaction(c,ca,cp,delta,p):
    """Equations 10-11 plus the accepted-increment Q defined in Section 2.2."""
    ra=np.where(ca<=p.transition_a,p.beta_ar*c-p.lambda_a*ca,
                p.beta_aa*c-(p.beta_ad+p.lambda_a)*ca)
    rp=np.where(cp<=p.transition_p,p.beta_p0*c-p.lambda_p*cp,
                p.beta_p0*p.transition_p*c/np.maximum(cp,p.transition_p)-p.lambda_p*cp)
    ra=np.where(ca>=p.capacity_a,np.minimum(ra,0),ra)
    rp=np.where(cp>=p.capacity_p,np.minimum(rp,0),rp)
    ta=ca+delta*ra;tp=cp+delta*rp
    an=np.clip(ta,0,p.capacity_a);pn=np.clip(tp,0,p.capacity_p)
    q=((an-ca)+(pn-cp))/delta+p.lambda_a*ca+p.lambda_p*cp
    clipping=int(np.count_nonzero(ta!=an)+np.count_nonzero(tp!=pn))
    return an,pn,q/p.porosity,clipping

def inlet_value(p,kind,t):
    if kind=="constant":return p.c_inlet
    if kind.startswith("exp_"):return p.c_inlet*(-math.expm1(-float(kind[4:])*t))
    raise ValueError(kind)

class Transport:
    def __init__(self,g,dispersion="constant",inlet="constant"):
        self.g=g;self.dispersion=dispersion;self.inlet=inlet;p=g.p
        self.c=np.zeros(g.shape);self.ca=np.zeros(g.shape);self.cp=np.zeros(g.shape)
        self.inlet_a=np.array(0.);self.inlet_p=np.array(0.)
        if dispersion=="constant":
            self.Dx=np.full((g.I,g.J-1),p.D_x)
            self.Dy=np.full((g.I-1,g.J),p.D_y)
        elif dispersion=="variable":
            xx,yy=np.meshgrid(g.x,g.y,indexing="ij")
            D=p.D0*(1+np.exp(-p.omega_D*xx))*(1+4*yy/p.H*(1-yy/p.H))
            self.Dx=(D[:-1,1:-1]+D[1:,1:-1])/2
            self.Dy=(D[1:-1,:-1]+D[1:-1,1:])/2
        else:raise ValueError(dispersion)
        self.cum_in=self.cum_out=self.cum_decay=0.
        self.max_half_residual=0.;self.max_cumulative_residual=0.
        self.min_c=0.;self.max_c=0.;self.clip_count=0
        self.max_symmetry=0.
    def inventory(self):
        p=self.g.p
        return self.g.cell_area*float(np.sum(p.porosity*self.c+self.ca+self.cp))
    def step(self,t,wx,wy0,wy1):
        p=self.g.p;delta=p.dt/2;theta=p.porosity;g=self.g
        X=Direction(g,0,wx/theta,self.Dx)
        Y0=Direction(g,1,wy0/theta,self.Dy)
        Y1=Direction(g,1,wy1/theta,self.Dy)
        c0=self.c; a0=self.ca;b0=self.cp
        inlet_half=inlet_value(p,self.inlet,t+delta)
        aa,bb,Q,clips=reaction(c0,a0,b0,delta,p)
        rhs=c0+delta*Y0.apply(c0)-delta*Q+delta*X.source(inlet_half)
        cs=X.solve(rhs,delta,p.lambda_c)
        fx=X.flux(cs,inlet_half);fy=Y0.flux(c0)
        self.record_half(c0,a0,b0,cs,aa,bb,fx,fy,clips)
        an,bn,Q,clips=reaction(cs,aa,bb,delta,p)
        cn=Y1.solve(cs+delta*X.apply(cs,inlet_half)-delta*Q,delta,p.lambda_c)
        fy=Y1.flux(cn)
        self.record_half(cs,aa,bb,cn,an,bn,fx,fy,clips)
        self.c=cn;self.ca=an;self.cp=bn
        for ct in (inlet_value(p,self.inlet,t),inlet_half):
            self.inlet_a,self.inlet_p,_,_=reaction(np.array(ct),self.inlet_a,self.inlet_p,delta,p)
        residual=self.inventory()-(self.cum_in-self.cum_out-self.cum_decay)
        self.max_cumulative_residual=max(self.max_cumulative_residual,abs(residual))
        self.max_symmetry=max(self.max_symmetry,float(np.max(np.abs(cn-cn[:,::-1]))))
        return residual
    def record_half(self,c0,a0,b0,cn,an,bn,fx,fy,clips):
        g=self.g;p=g.p;delta=p.dt/2;vol=g.cell_area;theta=p.porosity
        bi,bo=pore_boundary(fx,fy,g);bi*=theta*delta;bo*=theta*delta
        loss=delta*vol*float(np.sum(theta*p.lambda_c*cn+p.lambda_a*a0+p.lambda_p*b0))
        change=vol*float(np.sum(theta*(cn-c0)+(an-a0)+(bn-b0)))
        res=change-bi+bo+loss
        self.max_half_residual=max(self.max_half_residual,abs(res))
        self.cum_in+=bi;self.cum_out+=bo;self.cum_decay+=loss
        self.clip_count+=clips
        self.min_c=min(self.min_c,float(cn.min()))
        self.max_c=max(self.max_c,float(cn.max()))
        if not np.isfinite(cn).all() or cn.min() < -1e-12:
            raise ArithmeticError(f"Nonphysical dissolved solution: min={cn.min()}, dt={p.dt}")
    def snapshot(self,t):
        g=self.g
        return dict(c=g.full(self.c,inlet_value(g.p,self.inlet,t)),
                    ca=g.full(self.ca,float(self.inlet_a)),
                    cp=g.full(self.cp,float(self.inlet_p)))
    def metrics(self):
        return dict(inventory=self.inventory(),cumulative_in=self.cum_in,cumulative_out=self.cum_out,
                    cumulative_degradation=self.cum_decay,
                    final_residual=self.inventory()-self.cum_in+self.cum_out+self.cum_decay,
                    max_cumulative_residual=self.max_cumulative_residual,
                    relative_max_cumulative_residual=self.max_cumulative_residual/max(self.cum_in,1e-300),
                    max_half_step_residual=self.max_half_residual,min_c=self.min_c,max_c=self.max_c,
                    max_ca=float(self.ca.max()),max_cp=float(self.cp.max()),
                    capacity_projection_count=self.clip_count,max_symmetry_error=self.max_symmetry)

def simulate(p,output,scenarios=None,snapshot_times=(0,1800,3600,7200,18000,54000),log_every=3600.):
    """Save reproducible compact outputs; never loads any old result array."""
    output=Path(output)
    if output.exists():raise FileExistsError(f"Refusing to overwrite run {output}")
    output.mkdir(parents=True)
    g=Grid(p);flow=Pressure(g)
    if scenarios is None:scenarios=[("constant","constant"),("variable","constant")]
    transports={disp if inlet=="constant" else inlet:Transport(g,disp,inlet) for disp,inlet in scenarios}
    n=int(round(p.duration/p.dt));assert abs(n*p.dt-p.duration)<1e-8
    sample_stride=max(1,int(round(10/p.dt)))
    snapshots={key:{} for key in transports};flow_snaps={}
    events={int(round(t/p.dt)):t for t in snapshot_times if t<=p.duration}
    events[n]=p.duration
    history={key:[] for key in transports}
    flow_res=0.;min_p=p.p_ambient;max_p=p.p_ambient;flow_sym=0.
    started=time.perf_counter()
    def save_state(k,t):
        if k in events:
            pp=g.full(flow.q,flow.jump)+p.p_ambient
            wx,wy=flow.velocities(flow.q)
            flow_snaps[str(k)]=dict(t=np.array(t),p=pp,wx=wx,wy=wy)
            for key,tr in transports.items():snapshots[key][str(k)]=dict(t=np.array(t),**tr.snapshot(t))
        if k%sample_stride==0 or k==n:
            for key,tr in transports.items():
                m=tr.metrics()
                history[key].append([t,m['inventory'],m['cumulative_in'],m['cumulative_out'],
                    m['cumulative_degradation'],m['final_residual'],m['min_c'],
                    tr.c[0,g.J//2-1],tr.ca[0,g.J//2-1],tr.cp[0,g.J//2-1]])
    save_state(0,0)
    for k in range(n):
        wx,wy0,wy1,res=flow.step()
        flow_res=max(flow_res,abs(res))
        q=flow.q
        min_p=min(min_p,p.p_ambient+float(q.min()))
        max_p=max(max_p,p.p_ambient+float(q.max()))
        flow_sym=max(flow_sym,float(np.max(np.abs(q-q[:,::-1]))))
        for tr in transports.values():tr.step(k*p.dt,wx,wy0,wy1)
        save_state(k+1,(k+1)*p.dt)
        if (k+1)%max(1,int(round(log_every/p.dt)))==0:
            print(f"{output.name}: t={(k+1)*p.dt/3600:.2f} h, elapsed={time.perf_counter()-started:.1f} s",flush=True)
    np.savez_compressed(output/'flow_snapshots.npz',x=g.x,y=g.y,
        **{f'{key}_{field}':v for key,d in flow_snaps.items() for field,v in d.items()})
    for key,ss in snapshots.items():
        np.savez_compressed(output/(key+'_snapshots.npz'),x=g.x,y=g.y,
            **{f'{step}_{field}':v for step,d in ss.items() for field,v in d.items()})
        np.savetxt(output/(key+'_history.csv'),history[key],delimiter=',',
            header='time_s,inventory,cumulative_in,cumulative_out,cumulative_degradation,residual,min_c,c_at_dx_mid,ca_at_dx_mid,cp_at_dx_mid',comments='')
    summary={'parameters':asdict(p),'scenario_names':list(transports),'elapsed_seconds':time.perf_counter()-started,
        'steps':n,'snapshot_steps':events,'inventory_definition':'dx*dy*sum(theta*c+ca+cp) over interior solved control volumes; boundary nodes excluded',
        'inventory_area':(p.L-p.dx)*(p.H-p.dy),'effective_inlet_width':float(g.inlet_fraction.sum()*p.dy),
        'flow':{'max_half_pair_storage_residual':flow_res,'min_pressure':min_p,'max_interior_pressure':max_p,
                'max_symmetry_error':flow_sym,'storage_change_relative_to_reference_porosity':p.storage*(p.p_inlet-p.p_ambient)/p.porosity},
        'transport':{key:tr.metrics() for key,tr in transports.items()}}
    (output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary['transport'],indent=2),flush=True)
    return summary

