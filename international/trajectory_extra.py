"""Compute the unprotected full-score-space diagnostic after trajectory fitting."""
import argparse,sys
from pathlib import Path
from fig5_chain import MODELS,read_cortex,model_inputs,make,artifact_endpoint,fit
from trajectory_models import fingerprint
import numpy as np,json
from scipy.io import savemat


def extra(s,c):
    a,_=artifact_endpoint(s,c,s['pool'],'screen_r001');n=a['n'];ng=np.sqrt(n);Y=a['Y'].copy();ranks=[]
    for j in range(Y.shape[1]):
        Q,v,_=np.linalg.svd(ng[:,None]*np.c_[a['U'][j].T,np.ones(len(n)),c],full_matrices=False)
        Q=Q[:,v>v[0]*1e-12];rm=Q@(Q.T@(a['zy'][:,j]/ng));Y[:,j]-=(rm/ng)[a['g']];ranks.append(Q.shape[1])
    b,_=fit(make(a['X'],Y,dict(g=a['g'],C=c,targets=a['targets'])),c,0)
    norm=np.linalg.norm(b);valid=norm>1e-10
    r=np.array([np.corrcoef(b,v)[0,1] for v in a['targets'].T]) if valid else np.full(a['targets'].shape[1],np.nan)
    return b,r,norm,np.array(ranks),valid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['inputs','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--harmonized',type=Path);p.add_argument('--draws',type=int,default=0);p.add_argument('--seed',type=int,default=2026102201)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    for j,model in enumerate(MODELS):
        cohort=MODELS[model][0]
        if args.harmonized:
            with np.load(args.harmonized/(cohort+'_harmonized.npz')) as f:raw=dict(f)
        else:raw=read_cortex(args.inputs,cohort)
        d,*_=model_inputs(raw,args.inputs,model);s=make(d['X'],d['Y'],d);s['pool']=d['pool']
        rng=np.random.default_rng(args.seed+j);orders=np.array([rng.permutation(len(d['C'])) for _ in range(args.draws)])
        results=[extra(s,d['C'])]+[extra(s,d['C'][order]) for order in orders]
        b,r,norm,rank,valid=map(np.array,zip(*results))
        data=dict(beta=b,r=r,map_norm=norm,score_space_rank=rank,valid_spatial_correlation=valid,permutations=orders,roi_names=d['roi_names'],input_sha256=fingerprint(d),seed=args.seed+j)
        np.savez_compressed(args.output/(model+'_extra.npz'),**data);savemat(args.output/(model+'_extra.mat'),data,do_compression=True)
        print(model,'valid maps',int(valid.sum()),'/',len(valid),'observed norm',norm[0],flush=True)
