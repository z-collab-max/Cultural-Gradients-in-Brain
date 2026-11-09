"""Evaluate protected fits at every artifact count using cross-products."""
import numpy as np
from protected_reml import _legacy_fit


def fit_levels(s,c,levels=None):
    base,anchor=_legacy_fit(s,c,0)
    n=s['n'];N=s['N'];xz=s['xz'];xx=s['xx'];R=len(anchor)
    w=anchor[:,None]/(1+anchor[:,None]*n)
    vx=xx[None,:,:]-np.einsum('pg,rg,qg->rpq',xz,w,xz)
    vc=(xz@c)[None,:]-np.einsum('pg,rg,g->rp',xz,w,n*c)
    a=np.linalg.solve(vx,vc)
    t=c[None,:]-w*(n*c[None,:]-a@xz)
    zh=t*n-a@xz
    xh=t@xz.T-a@xx
    hh=np.sum(t*t*n,axis=1)-2*np.sum((a@xz)*t,axis=1)+np.einsum('rp,pq,rq->r',a,xx,a)
    if np.any(hh<=0):raise ValueError('Degenerate protected cultural direction')
    hn=np.sqrt(hh)
    zq=np.stack([np.broadcast_to(n/np.sqrt(N),(R,len(n))),zh/hn[:,None]],axis=2)
    xq=np.stack([np.broadcast_to(s['X'].sum(0)/np.sqrt(N),(R,xx.shape[0])),xh/hn[:,None]],axis=2)
    qy=np.column_stack([s['Y'].sum(0)/np.sqrt(N),(np.sum(t*s['zy'].T,axis=1)-np.sum(a*s['xy'].T,axis=1))/hn])
    if levels is None:levels=range(s['U'].shape[1]+1)
    maps=[]
    for k in levels:
        if k==0:
            maps.append(base);continue
        U=s['U'][:,:k];h=U@zq
        rhs=np.einsum('rkg,gr->rk',U,s['zy'])-np.einsum('rkl,rl->rk',h,qy)
        gram=np.einsum('rkg,g,rlg->rkl',U,n,U)-h@h.transpose(0,2,1)
        ev,V=np.linalg.eigh((gram+gram.transpose(0,2,1))/2)
        inv=np.zeros_like(ev);np.divide(1,ev,out=inv,where=ev>1e-10)
        coef=np.einsum('rki,ri->rk',V,np.einsum('rki,rk->ri',V,rhs)*inv)
        rm=np.einsum('rkg,rk->gr',U,coef);qc=np.einsum('rkl,rk->rl',h,coef)
        q=dict(s,xy=s['xy']-xz@rm+np.einsum('rpl,rl->pr',xq,qc),
               zy=s['zy']-n[:,None]*rm+np.einsum('rgl,rl->gr',zq,qc),
               yy=s['yy']-np.sum(rhs*coef,axis=1))
        maps.append(_legacy_fit(q,c,0)[0])
    return np.array(maps)
