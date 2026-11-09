"""Raw cortical chain to all-PC trajectories and unrestricted shuffle tests."""
import argparse,hashlib,json
from pathlib import Path
import sys
from fig5_chain import MODELS,read_cortex,model_inputs,make,artifact_endpoint,fit,unit
from trajectory_fit import fit_levels
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
from panel_data import panel,save


def fingerprint(d):
    h=hashlib.sha256()
    for key in ['X','Y','g','C','pool','targets']:
        a=np.ascontiguousarray(d[key]);h.update(key.encode());h.update(str((a.shape,str(a.dtype))).encode());h.update(a.tobytes())
    return h.hexdigest()


def run(d,model,out,draws,seed,verify):
    cache=out/model/'candidates';cache.mkdir(parents=True,exist_ok=True)
    signature=dict(input_sha256=fingerprint(d),seed=seed,draws=draws,algorithm='protected_crossproducts_v1')
    path=cache/'configuration.json'
    if path.exists() and json.loads(path.read_text())!=signature:raise ValueError('Choose a new output directory for different inputs or settings')
    path.write_text(json.dumps(signature,indent=2))
    s=make(d['X'],d['Y'],d);rng=np.random.default_rng(seed)
    orders=np.array([rng.permutation(len(d['C'])) for _ in range(draws)])
    values=[];maps0=None
    for j in range(draws+1):
        f=cache/f'candidate_{j:05d}.npz'
        if f.exists():
            with np.load(f) as old:beta=old['beta'];r=old['r']
        else:
            c=d['C'] if j==0 else d['C'][orders[j-1]]
            a,retained=artifact_endpoint(s,c,d['pool'],'screen_r001');beta=fit_levels(a,c)
            if verify and j==0:
                for k in [0,1,a['k']//2,a['k']]:
                    reference=fit(a,c,k)[0]
                    np.testing.assert_allclose(beta[k],reference,atol=1e-7,rtol=1e-6)
            r=unit(beta.T).T@unit(d['targets'])
            np.savez_compressed(f,beta=beta,r=r,artifact_pool_count=retained)
        if j==0:maps0=beta
        values.append(r)
        if j%10==0 or j==draws:print(model,'trajectory',j,'/',draws,flush=True)
    values=np.array(values);obs=values[0];null=values[1:]
    p=(abs(null)>=abs(obs)).mean(0)
    source=dict(r=obs,null_r=null,p=p,observed_maps=maps0,permutations=orders,roi_names=d['roi_names'],removed_pcs=np.arange(len(obs)),endpoint_k=len(obs)-1)
    panels=[]
    for t,name in enumerate(MODELS[model][1:3]):
        panels.append(panel('ExtendedDataFigure2',model+'_'+name+'_trajectory','trajectory',r=obs[:,t],p=p[:,t],selected=len(obs)-1,**{f'null_{k}':null[:,k,t] for k in range(len(obs))}))
    save('ExtendedDataFigure2',out/model,source,panels)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['inputs','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--harmonized',type=Path);p.add_argument('--draws',type=int,default=1000)
    p.add_argument('--seed',type=int,default=2026102201);p.add_argument('--models',nargs='+',choices=MODELS,default=list(MODELS));p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    if a.draws<2:p.error('At least two shuffles are required')
    for model in a.models:
        cohort=MODELS[model][0]
        if a.harmonized:
            with np.load(a.harmonized/(cohort+'_harmonized.npz')) as f:raw=dict(f)
        else:raw=read_cortex(a.inputs,cohort)
        d,*_=model_inputs(raw,a.inputs,model);run(d,model,a.output,a.draws,a.seed+list(MODELS).index(model),a.verify)
