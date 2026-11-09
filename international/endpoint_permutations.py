"""Unit shuffles with endpoint refits and restartable block output."""
import hashlib,json
from pathlib import Path
import numpy as np
from artifact import artifact_endpoint
from protected_reml import fit
from trajectory_fit import fit_levels
from core import unit


def permutation_test(state,d,observed,repetitions,seed,cache=None,*,orders=None):
    supplied=orders is not None
    if supplied:
        orders=np.asarray(orders)
        if orders.dtype.kind not in 'iu' or orders.shape!=(repetitions,len(d['C'])) or not np.all(np.sort(orders,axis=1)==np.arange(len(d['C']))):
            raise ValueError('Expected zero-based draw-by-unit permutation indices')
        orders=orders.astype(np.int64,copy=True)
    else:
        rng=np.random.default_rng(seed)
        orders=np.array([rng.permutation(len(d['C'])) for _ in range(repetitions)])
    method='unrestricted unit shuffle'
    if cache is not None:
        cache=Path(cache);cache.mkdir(parents=True,exist_ok=True);h=hashlib.sha256()
        for key in ['X','Y','g','C','pool','targets']:
            a=np.ascontiguousarray(d[key]);h.update(key.encode());h.update(str((a.shape,str(a.dtype))).encode());h.update(a.tobytes())
        signature=dict(inputs=h.hexdigest(),draws=repetitions,seed=None if supplied else seed,
                       algorithm='endpoint_direct_absolute_v2')
        if supplied:signature['orders_sha256']=hashlib.sha256(np.ascontiguousarray(orders).tobytes()).hexdigest()
        config=cache/'configuration.json'
        if config.exists() and json.loads(config.read_text())!=signature:raise ValueError('Cached permutation settings differ')
        config.write_text(json.dumps(signature,indent=2))
    maps=[];counts=[]
    for first in range(0,repetitions,100):
        path=cache/f'block_{first:05d}.npz' if cache is not None else None
        if path is not None and path.exists():
            with np.load(path) as f:beta=f['beta'];retained=f['retained']
        else:
            beta=[];retained=[]
            for j,order in enumerate(orders[first:first+100]):
                c=d['C'][order];a,n=artifact_endpoint(state,c,d['pool'],'screen_r001')
                b=fit_levels(a,c,[a['k']])[0]
                if first+j<3:np.testing.assert_allclose(b,fit(a,c,a['k'])[0],atol=1e-7,rtol=1e-6)
                beta.append(b);retained.append(n)
            beta=np.array(beta);retained=np.array(retained)
            if path is not None:
                tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,beta=beta,retained=retained);tmp.replace(path)
        maps.append(beta);counts.append(retained)
        if first%1000==0:print('Endpoint permutations',first+len(beta),'/',repetitions,flush=True)
    maps=np.concatenate(maps);null_r=unit(maps.T).T@unit(d['targets'])
    p=(abs(null_r)>=abs(observed)).mean(0)
    return dict(permutations=orders,null_beta=maps,null_r=null_r,p_two_sided_absolute=p,
                artifact_pool_counts=np.concatenate(counts),null_method=method)
