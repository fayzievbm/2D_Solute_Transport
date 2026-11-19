"""Compiled conservative ADI kernels, float64; fastmath is deliberately disabled."""
import numpy as np
from numba import njit

@njit(cache=True)
def coefficients(v,D,h,axis,fraction):
    a=np.empty_like(v);b=np.empty_like(v)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            d=D[i,j]/(h[i] if axis==0 else h[j])
            if axis==0 and i==0:d*=fraction[j]
            if axis==1 and (j==0 or j==v.shape[1]-1):d=0.
            a[i,j]=max(v[i,j],0.)+d
            b[i,j]=min(v[i,j],0.)-d
    return a,b

@njit(cache=True)
def flux(z,a,b,axis,left=0.):
    out=np.zeros_like(a);nx,ny=z.shape
    if axis==0:
        for j in range(ny):
            out[0,j]=a[0,j]*left+b[0,j]*z[0,j]
            out[nx,j]=a[nx,j]*z[nx-1,j]
        for i in range(1,nx):
            for j in range(ny):out[i,j]=a[i,j]*z[i-1,j]+b[i,j]*z[i,j]
    else:
        for i in range(nx):
            for j in range(1,ny):out[i,j]=a[i,j]*z[i,j-1]+b[i,j]*z[i,j]
    return out

@njit(cache=True)
def apply(z,a,b,axis,width,left=0.):
    ff=flux(z,a,b,axis,left);out=np.empty_like(z)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            out[i,j]=(-(ff[i+1,j]-ff[i,j])/width[i] if axis==0 else
                      -(ff[i,j+1]-ff[i,j])/width[j])
    return out

@njit(cache=True)
def source(z,a,delta,width,left):
    rhs=z.copy()
    for j in range(z.shape[1]):rhs[0,j]+=delta*a[0,j]*left/width[0]
    return rhs

@njit(cache=True)
def solve(rhs,a,b,axis,width,delta,decay=0.):
    nx,ny=rhs.shape;out=rhs.copy();up=np.empty_like(rhs)
    if axis==0:
        for i in range(nx):
            for j in range(ny):
                lo=-delta*a[i,j]/width[i]
                diag=1+delta*decay-delta*(b[i,j]-a[i+1,j])/width[i]
                uu=delta*b[i+1,j]/width[i]
                if i>0:
                    diag-=lo*up[i-1,j];out[i,j]-=lo*out[i-1,j]
                up[i,j]=uu/diag;out[i,j]/=diag
        for i in range(nx-2,-1,-1):
            for j in range(ny):out[i,j]-=up[i,j]*out[i+1,j]
    else:
        for i in range(nx):
            for j in range(ny):
                lo=-delta*a[i,j]/width[j]
                diag=1+delta*decay-delta*(b[i,j]-a[i,j+1])/width[j]
                uu=delta*b[i,j+1]/width[j]
                if j>0:
                    diag-=lo*up[i,j-1];out[i,j]-=lo*out[i,j-1]
                up[i,j]=uu/diag;out[i,j]/=diag
            for j in range(ny-2,-1,-1):out[i,j]-=up[i,j]*out[i,j+1]
    return out

@njit(cache=True)
def react(c,a,b,delta,r):
    # r: theta, c0, ca0, cp0, car, cp1, bar, baa, bad, bp, lc, la, lp
    theta,c0,ca0,cp0,car,cp1,bar,baa,bad,bp,lc,la,lp=r
    ra=(bar*c-la*a) if a<=car else (baa*c-(bad+la)*a)
    rp=(bp*c-lp*b) if b<=cp1 else (bp*cp1*c/b-lp*b)
    if a>=ca0:ra=min(ra,0.)
    if b>=cp0:rp=min(rp,0.)
    ta=a+delta*ra;tb=b+delta*rp
    an=min(max(ta,0.),ca0);bn=min(max(tb,0.),cp0)
    return an,bn,((an-a)+(bn-b))/delta/theta+(la*a+lp*b)/theta,int(ta!=an)+int(tb!=bn)

@njit(cache=True)
def reactions(c,a,b,delta,r):
    an=np.empty_like(c);bn=np.empty_like(c);Q=np.empty_like(c);clips=0
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):
            an[i,j],bn[i,j],Q[i,j],n=react(c[i,j],a[i,j],b[i,j],delta,r);clips+=n
    return an,bn,Q,clips

@njit(cache=True)
def amount(c,a,b,volx,voly,theta):
    m=0.
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):m+=volx[i]*voly[j]*(theta*c[i,j]+a[i,j]+b[i,j])
    return m

@njit(cache=True)
def ledger(c,a,b,cn,an,bn,fx,fy,volx,voly,delta,r,stats,clips):
    bi=0.;bo=0.;loss=0.;change=0.;theta=r[0]
    for j in range(fx.shape[1]):
        left=fx[0,j]*voly[j]*theta*delta;right=-fx[-1,j]*voly[j]*theta*delta
        bi+=max(left,0)+max(right,0);bo-=min(left,0)+min(right,0)
    for i in range(fy.shape[0]):
        low=fy[i,0]*volx[i]*theta*delta;high=-fy[i,-1]*volx[i]*theta*delta
        bi+=max(low,0)+max(high,0);bo-=min(low,0)+min(high,0)
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):
            v=volx[i]*voly[j]
            loss+=delta*v*(theta*r[10]*cn[i,j]+r[11]*a[i,j]+r[12]*b[i,j])
            change+=v*(theta*(cn[i,j]-c[i,j])+an[i,j]-a[i,j]+bn[i,j]-b[i,j])
            stats[5]=min(stats[5],cn[i,j]);stats[6]=max(stats[6],cn[i,j])
            if not np.isfinite(cn[i,j]) or cn[i,j]<-1e-12:raise ValueError('Negative/nonfinite concentration')
    stats[0]+=bi;stats[1]+=bo;stats[2]+=loss
    stats[4]=max(stats[4],abs(change-bi+bo+loss));stats[7]+=clips
    residual=amount(cn,an,bn,volx,voly,theta)-stats[0]+stats[1]+stats[2]
    stats[3]=max(stats[3],abs(residual))

@njit(cache=True)
def advance(q,states,reservoir,stats,flowstats,kstart,nsteps,dt,hx,hy,volx,voly,
            fraction,pxa,pxb,pya,pyb,Dsx,Dsy,phis,r,jump,storage):
    delta=dt/2
    for k in range(kstart,kstart+nsteps):
        rhs=source(q+delta*apply(q,pya,pyb,1,voly),pxa,delta,volx,jump)
        qs=solve(rhs,pxa,pxb,0,volx,delta)
        qn=solve(qs+delta*apply(qs,pxa,pxb,0,volx,jump),pya,pyb,1,voly,delta)
        wx=storage*flux(qs,pxa,pxb,0,jump)
        wy0=storage*flux(q,pya,pyb,1);wy1=storage*flux(qn,pya,pyb,1)
        vchange=0.;net=0.
        for i in range(q.shape[0]):
            for j in range(q.shape[1]):
                vchange+=storage*volx[i]*voly[j]*(qn[i,j]-q[i,j])
                flowstats[1]=min(flowstats[1],qn[i,j]);flowstats[2]=max(flowstats[2],qn[i,j])
        for j in range(wx.shape[1]):net+=dt*(wx[0,j]-wx[-1,j])*voly[j]
        flowstats[0]=max(flowstats[0],abs(vchange-net));q=qn
        for s in range(states.shape[0]):
            c=states[s,0];a=states[s,1];b=states[s,2]
            xa,xb=coefficients(wx/r[0],Dsx[s],hx,0,fraction)
            ya0,yb0=coefficients(wy0/r[0],Dsy[s],hy,1,fraction)
            ya1,yb1=coefficients(wy1/r[0],Dsy[s],hy,1,fraction)
            ih=r[1] if phis[s]==0 else r[1]*(-np.expm1(-phis[s]*(k*dt+delta)))
            aa,bb,Q,clips=reactions(c,a,b,delta,r)
            rhs=source(c+delta*apply(c,ya0,yb0,1,voly)-delta*Q,xa,delta,volx,ih)
            cs=solve(rhs,xa,xb,0,volx,delta,r[10])
            fx=flux(cs,xa,xb,0,ih);fy=flux(c,ya0,yb0,1)
            ledger(c,a,b,cs,aa,bb,fx,fy,volx,voly,delta,r,stats[s],clips)
            an,bn,Q,clips=reactions(cs,aa,bb,delta,r)
            cn=solve(cs+delta*apply(cs,xa,xb,0,volx,ih)-delta*Q,ya1,yb1,1,voly,delta,r[10])
            fy=flux(cn,ya1,yb1,1)
            ledger(cs,aa,bb,cn,an,bn,fx,fy,volx,voly,delta,r,stats[s],clips)
            states[s,0]=cn;states[s,1]=an;states[s,2]=bn
            it=r[1] if phis[s]==0 else r[1]*(-np.expm1(-phis[s]*k*dt))
            for val in (it,ih):
                reservoir[s,0],reservoir[s,1],unused,unused2=react(val,reservoir[s,0],reservoir[s,1],delta,r)
    return q
