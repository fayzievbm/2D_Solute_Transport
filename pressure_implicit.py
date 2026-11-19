"""Unsplit finite-volume pressure ODE, implicit BDF, exact structured solves.

The finite inlet is a diagonal perturbation on the first interior x row.
A transverse symmetric eigendecomposition and a Woodbury correction solve
the FULL implicit matrix; this is not directional time splitting.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.linalg import eigh_tridiagonal,lu_factor,lu_solve
from scipy.integrate import BDF
from numba import njit

@njit(cache=True)
def mode_factors(lower,diagonal,upper,eigenvalues,gamma):
    nx=len(diagonal);ny=len(eigenvalues)
    inv=np.empty((nx,ny));up=np.empty((nx,ny));green=np.empty((nx,ny))
    for i in range(nx):
        lo=-gamma*lower[i];uu=-gamma*upper[i]
        for j in range(ny):
            dd=1-gamma*(diagonal[i]+eigenvalues[j])
            if i:dd-=lo*up[i-1,j]
            inv[i,j]=1/dd;up[i,j]=uu/dd
            green[i,j]=(1. if i==0 else -lo*green[i-1,j])*inv[i,j]
    for i in range(nx-2,-1,-1):
        for j in range(ny):green[i,j]-=up[i,j]*green[i+1,j]
    return inv,up,green

@njit(cache=True)
def mode_solve(rhs,lower,gamma,inv,up):
    z=rhs.copy();nx,ny=z.shape
    for i in range(nx):
        for j in range(ny):
            if i:z[i,j]+=gamma*lower[i]*z[i-1,j]
            z[i,j]*=inv[i,j]
    for i in range(nx-2,-1,-1):
        for j in range(ny):z[i,j]-=up[i,j]*z[i+1,j]
    return z

class PressureOperator:
    def __init__(self,m):
        p=m.p;self.p=p;self.nx=len(m.vx);self.ny=len(m.vy)
        self.shape=(self.nx,self.ny);self.vx=m.vx;self.vy=m.vy
        self.volume=(m.vx[:,None]*m.vy[None,:]).ravel()
        ax=p.permeability_x/(p.viscosity*p.storage)
        ay=p.permeability_y/(p.viscosity*p.storage)
        gx=ax/m.hx;gy=ay/m.hy;gy[0]=gy[-1]=0.
        self.lower=gx[:-1]/m.vx;self.lower[0]=0.
        self.upper=gx[1:]/m.vx;self.diagonal=-(self.lower+self.upper)
        self.upper[-1]=0. # right boundary remains in diagonal
        ydiag=-(gy[:-1]+gy[1:])/m.vy
        yoff=gy[1:-1]/np.sqrt(m.vy[:-1]*m.vy[1:])
        self.eig,self.Q=eigh_tridiagonal(ydiag,yoff,lapack_driver='stev')
        self.sqrtvy=np.sqrt(m.vy)
        self.beta=gx[0]/m.vx[0]*m.f
        self.open=np.flatnonzero(self.beta>0)
        self.Qopen=self.Q[self.open,:]
        ids=np.arange(self.nx*self.ny).reshape(self.shape)
        diag=self.diagonal[:,None]+ydiag[None,:]
        diag=diag.copy();diag[0]-=self.beta
        rows=np.r_[ids.ravel(),ids[1:].ravel(),ids[:-1].ravel(),ids[:,1:].ravel(),ids[:,:-1].ravel()]
        cols=np.r_[ids.ravel(),ids[:-1].ravel(),ids[1:].ravel(),ids[:,:-1].ravel(),ids[:,1:].ravel()]
        values=np.r_[diag.ravel(),np.broadcast_to(self.lower[1:,None],ids[1:].shape).ravel(),
            np.broadcast_to(self.upper[:-1,None],ids[:-1].shape).ravel(),
            np.broadcast_to((gy[1:-1]/m.vy[1:])[None,:],ids[:,1:].shape).ravel(),
            np.broadcast_to((gy[1:-1]/m.vy[:-1])[None,:],ids[:,:-1].shape).ravel()]
        self.A=coo_matrix((values,(rows,cols)),shape=(ids.size,ids.size)).tocsc()
        self.Acsr=self.A.tocsr()
        self.b=np.zeros(self.shape);self.b[0]=self.beta*(p.p_inlet-p.p_ambient);self.b=self.b.ravel()
        self.flow_left=p.storage*gx[0]*m.f*m.vy
        self.flow_right=-p.storage*gx[-1]*m.vy
        self.flow_constant=float(self.flow_left.sum()*(p.p_inlet-p.p_ambient))
        # Warm both compiled signatures before setting hardware FP mode.
        f=self.factor(1.);self.solve(f,np.zeros(ids.size))

    def factor(self,gamma):
        gamma=float(gamma)
        inv,up,H=mode_factors(self.lower,self.diagonal,self.upper,self.eig,gamma)
        D=gamma*self.beta[self.open]
        G=(self.Qopen*H[0])@self.Qopen.T
        small=np.eye(len(self.open))+G*D[None,:]
        return gamma,inv,up,H,D,lu_factor(small,check_finite=False)

    def solve(self,factor,rhs):
        gamma,inv,up,H,D,small=factor
        modes=(rhs.reshape(self.shape)*self.sqrtvy)@self.Q
        modes=mode_solve(modes,self.lower,gamma,inv,up)
        boundary=modes[0]@self.Qopen.T
        u=lu_solve(small,boundary,check_finite=False)
        modes-=H*(self.Qopen.T@(D*u))[None,:]
        return ((modes@self.Q.T)/self.sqrtvy).ravel()

    def rhs(self,t,q):return self.Acsr@q+self.b

    def flow_coefficients(self,dense):
        D=dense.D.reshape((-1,*self.shape))
        coefficients=-D[:,0,:]@self.flow_left+D[:,-1,:]@self.flow_right
        coefficients[0]+=self.flow_constant
        return coefficients

class Pressure:
    def __init__(self,m,duration,rtol=1e-7,atol=1e-7,linear_solver='structured'):
        self.op=PressureOperator(m);o=self.op
        self.solver=BDF(o.rhs,0.,np.zeros(o.A.shape[0]),float(duration),rtol=rtol,atol=atol,jac=o.A)
        if linear_solver=='structured':
            def factor(matrix):
                self.solver.nlu+=1
                gamma=-float(matrix[0,o.ny])/float(o.A[0,o.ny])
                return o.factor(gamma)
            self.solver.lu=factor;self.solver.solve_lu=o.solve
        self.dense=None;self.integral=0.;self.window_integral=0.;self.last_requested=0.
        self.stats={'min_gauge_pressure':0.,'max_gauge_pressure':0.,'adaptive_steps':0,
            'max_storage_integral_residual':0.,'max_relative_storage_integral_residual_after_1s':0.}
        self.rtol=rtol;self.atol=atol;self.linear_solver=linear_solver
        self.gauss_nodes,self.gauss_weights=np.polynomial.legendre.leggauss(3)

    def integrated_flux(self,t):
        left=self.dense.t_old;half=(t-left)/2
        times=left+half*(self.gauss_nodes+1)
        basis=np.vstack([np.ones(3),np.cumprod((times[None,:]-self.dense.t_shift[:,None])/self.dense.denom[:,None],axis=0)])
        return float(half*np.dot(self.gauss_weights,self.flux_coeff@basis))

    def check_bounds(self,q):
        lo=float(q.min());hi=float(q.max());jump=self.op.p.p_inlet-self.op.p.p_ambient
        self.stats['min_gauge_pressure']=min(self.stats['min_gauge_pressure'],lo)
        self.stats['max_gauge_pressure']=max(self.stats['max_gauge_pressure'],hi)
        if lo < -1e-6 or hi>jump+1e-6 or not np.isfinite(q).all():
            raise ValueError(f'Pressure bound failed: gauge range {lo:g} to {hi:g} Pa')

    def at(self,t):
        if t<self.last_requested-1e-9:raise ValueError('Pressure queries must be chronological')
        self.last_requested=t
        if t==0:return np.zeros(self.op.shape)
        while self.solver.t<t:
            if self.solver.status!='running':raise RuntimeError('Pressure integrator cannot advance')
            self.window_integral=self.integral
            message=self.solver.step()
            if self.solver.status=='failed':raise RuntimeError(message)
            self.dense=self.solver.dense_output();self.flux_coeff=self.op.flow_coefficients(self.dense)
            self.integral=self.window_integral+self.integrated_flux(self.solver.t)
            self.stats['adaptive_steps']+=1;self.check_bounds(self.solver.y)
        q=self.dense(t);self.check_bounds(q)
        integrated=self.window_integral+self.integrated_flux(t)
        storage=self.op.p.storage*float(np.dot(self.op.volume,q))
        residual=abs(storage-integrated)
        self.stats['max_storage_integral_residual']=max(self.stats['max_storage_integral_residual'],residual)
        if t>=1:self.stats['max_relative_storage_integral_residual_after_1s']=max(self.stats['max_relative_storage_integral_residual_after_1s'],residual/max(abs(integrated),1e-300))
        return q.reshape(self.op.shape)
