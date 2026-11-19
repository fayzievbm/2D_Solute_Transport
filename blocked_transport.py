"""Allocation-efficient equivalent ADI transport; no fastmath, no model change."""
from blocked_base import coefficients,flux,apply,source,solve,react,reactions,amount
import numpy as np
from numba import njit, prange

@njit(cache=True,inline='always')
def xc(i,j,w,D,h,f,theta):
    u=w[i,j]/theta;d=D[i,j]/h[i]
    if i==0:d*=f[j]
    return max(u,0.)+d,min(u,0.)-d

@njit(cache=True,inline='always')
def yc(i,j,w,D,h,theta):
    if j==0 or j==w.shape[1]-1:return 0.,0.
    u=w[i,j]/theta;d=D[i,j]/h[j]
    return max(u,0.)+d,min(u,0.)-d

@njit(cache=True, parallel=True)
def transport(c,a,b,wx,wy0,wy1,Dx,Dy,hx,hy,vx,vy,f,delta,ih,r,stats):
    nx,ny=c.shape;theta=r[0];lc=r[10];cs=np.empty_like(c);up=np.empty_like(c)
    change=0.;loss=0.;inventory=0.;clips=0;minc=stats[5];maxc=stats[6];invalid=0
    # X-implicit/Y-explicit stage. Reaction and source assembly are fused
    # with the Thomas forward sweep; mass accounting with backward sweep.
    for block in prange((ny+15)//16):
        jstart=block*16;jend=min(jstart+16,ny)
        for i in range(nx):
            for j in range(jstart,jend):
                volume=vx[i]*vy[j];co=c[i,j];ao=a[i,j];bo=b[i,j]
                an,bn,Q,nclip=react(co,ao,bo,delta,r);a[i,j]=an;b[i,j]=bn;clips+=nclip
                change+=volume*(an-ao+bn-bo);loss+=delta*volume*(r[11]*ao+r[12]*bo)
                al,bl=xc(i,j,wx,Dx,hx,f,theta);ar,br=xc(i+1,j,wx,Dx,hx,f,theta)
                ab,bb=yc(i,j,wy0,Dy,hy,theta);at,bt=yc(i,j+1,wy0,Dy,hy,theta)
                bottom=0. if j==0 else ab*c[i,j-1]+bb*co
                top=0. if j==ny-1 else at*co+bt*c[i,j+1]
                rhs=co-delta*(top-bottom)/vy[j]-delta*Q
                if i==0:rhs+=delta*al*ih/vx[i]
                lo=-delta*al/vx[i];diag=1+delta*lc-delta*(bl-ar)/vx[i];uu=delta*br/vx[i]
                if i>0:diag-=lo*up[i-1,j];rhs-=lo*cs[i-1,j]
                cs[i,j]=rhs/diag;up[i,j]=uu/diag
        for i in range(nx-1,-1,-1):
            for j in range(jstart,jend):
                if i<nx-1:cs[i,j]-=up[i,j]*cs[i+1,j]
                cn=cs[i,j];volume=vx[i]*vy[j]
                change+=volume*theta*(cn-c[i,j]);loss+=delta*volume*theta*lc*cn
                inventory+=volume*(theta*cn+a[i,j]+b[i,j])
                minc=min(minc,cn);maxc=max(maxc,cn)
                if not np.isfinite(cn) or cn<-1e-12:invalid+=1
    bi=0.;bout=0.
    for j in range(ny):
        al,bl=xc(0,j,wx,Dx,hx,f,theta);ar,br=xc(nx,j,wx,Dx,hx,f,theta)
        incoming=delta*theta*vy[j]*(al*ih+bl*cs[0,j]);outgoing=-delta*theta*vy[j]*ar*cs[-1,j]
        bi+=max(incoming,0)+max(outgoing,0);bout-=min(incoming,0)+min(outgoing,0)
    if invalid:raise ValueError('Negative/nonfinite concentration')
    stats[5]=minc;stats[6]=maxc
    stats[0]+=bi;stats[1]+=bout;stats[2]+=loss;stats[7]+=clips
    stats[3]=max(stats[3],abs(inventory-stats[0]+stats[1]+stats[2]))
    stats[4]=max(stats[4],abs(change-bi+bout+loss))
    change=0.;loss=0.;inventory=0.;clips=0;minc=stats[5];maxc=stats[6];invalid=0
    # Y-implicit/X-explicit stage. cs is immutable during this stage; c is output.
    for i in prange(nx):
        for j in range(ny):
            volume=vx[i]*vy[j];co=cs[i,j];ao=a[i,j];bo=b[i,j]
            an,bn,Q,nclip=react(co,ao,bo,delta,r);a[i,j]=an;b[i,j]=bn;clips+=nclip
            change+=volume*(an-ao+bn-bo);loss+=delta*volume*(r[11]*ao+r[12]*bo)
            al,bl=xc(i,j,wx,Dx,hx,f,theta);ar,br=xc(i+1,j,wx,Dx,hx,f,theta)
            left=al*(ih if i==0 else cs[i-1,j])+bl*co
            right=ar*co+(0. if i==nx-1 else br*cs[i+1,j])
            rhs=co-delta*(right-left)/vx[i]-delta*Q
            ab,bb=yc(i,j,wy1,Dy,hy,theta);at,bt=yc(i,j+1,wy1,Dy,hy,theta)
            lo=-delta*ab/vy[j];diag=1+delta*lc-delta*(bb-at)/vy[j];uu=delta*bt/vy[j]
            if j>0:diag-=lo*up[i,j-1];rhs-=lo*c[i,j-1]
            c[i,j]=rhs/diag;up[i,j]=uu/diag
        for j in range(ny-1,-1,-1):
            if j<ny-1:c[i,j]-=up[i,j]*c[i,j+1]
            cn=c[i,j];volume=vx[i]*vy[j]
            change+=volume*theta*(cn-cs[i,j]);loss+=delta*volume*theta*lc*cn
            inventory+=volume*(theta*cn+a[i,j]+b[i,j])
            minc=min(minc,cn);maxc=max(maxc,cn)
            if not np.isfinite(cn) or cn<-1e-12:invalid+=1
    if invalid:raise ValueError('Negative/nonfinite concentration')
    stats[5]=minc;stats[6]=maxc
    stats[0]+=bi;stats[1]+=bout;stats[2]+=loss;stats[7]+=clips
    stats[3]=max(stats[3],abs(inventory-stats[0]+stats[1]+stats[2]))
    stats[4]=max(stats[4],abs(change-bi+bout+loss))

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
            ih=r[1] if phis[s]==0 else r[1]*(-np.expm1(-phis[s]*(k*dt+delta)))
            transport(states[s,0],states[s,1],states[s,2],wx,wy0,wy1,Dsx[s],Dsy[s],hx,hy,volx,voly,fraction,delta,ih,r,stats[s])
            it=r[1] if phis[s]==0 else r[1]*(-np.expm1(-phis[s]*k*dt))
            for val in (it,ih):
                reservoir[s,0],reservoir[s,1],unused,unused2=react(val,reservoir[s,0],reservoir[s,1],delta,r)
    return q
