from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'core'))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analyses'))
import numpy as np
from core import z
from engine import make
from artifact import artifact_endpoint
from protected_reml import fit,_legacy_fit
from webb_statistics import spatial_statistics


def webb_intervals(d,repetitions=1000,seed=2026102201):
    if repetitions<2:raise ValueError('At least two Webb draws are required')
    for key in ['X','Y','C','targets']:
        if not np.isfinite(d[key]).all():raise ValueError('Nonfinite '+key)
    if not np.array_equal(np.unique(d['g']),np.arange(len(d['C']))):raise ValueError('g must contain contiguous zero-based group IDs')
    X,Y,g,c=d['X'],d['Y'],d['g'],d['C'];n=np.bincount(g);Z=np.eye(len(c))[g]
    design=np.c_[X,c[g]];mx=Z.T@design/n[:,None]
    state=make(X,Y,d);artifact,_=artifact_endpoint(state,c,d['pool'],'screen_r001')
    point,_,_=fit(artifact,c,artifact['k'],True)
    _,anchor=_legacy_fit(state,c,0);mu=np.zeros_like(Y);adjusted=np.zeros_like(Y)
    for j,lam in enumerate(anchor):
        alpha=1-1/np.sqrt(1+n*lam);xt=design-alpha[g,None]*mx[g]
        my=Z.T@Y[:,j]/n;yt=Y[:,j]-alpha[g]*my[g];Q,R=np.linalg.qr(xt,mode='reduced')
        b=np.linalg.solve(R,Q.T@yt);mu[:,j]=design@b;res=yt-xt@b;h=np.sum(Q*Q,axis=1)
        if np.max(h)>=1:raise ValueError('Unit leverage in bootstrap design')
        e=res/np.sqrt(1-h);inflate=np.sqrt(1+n*lam)-1
        adjusted[:,j]=e+inflate[g]*(Z.T@e/n)[g]
        np.testing.assert_allclose(mu[:,j]+res+inflate[g]*(Z.T@res/n)[g],Y[:,j],atol=1e-8)
    support=np.array([-np.sqrt(1.5),-1,-np.sqrt(.5),np.sqrt(.5),1,np.sqrt(1.5)])
    rng=np.random.RandomState(seed);weights=support[rng.randint(0,6,(repetitions,len(c)))]
    maps=[]
    for index,w in enumerate(weights,1):
        y=z(mu+adjusted*w[g,None]);s=make(X,y,d);a,_=artifact_endpoint(s,c,d['pool'],'screen_r001')
        beta,_=fit(a,c,a['k']);maps.append(beta)
        if index%100==0 or index==repetitions:print('Webb refits',index,'/',repetitions,flush=True)
    maps=np.asarray(maps);out=spatial_statistics(point,maps,d['targets'],quantile_method='linear')
    out.update(map_draws=maps,weights=weights,seed=seed)
    return out
