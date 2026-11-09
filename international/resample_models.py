"""Refit endpoint cortical maps after sampling individuals or cultural units."""
import argparse
import warnings
from pathlib import Path
import sys
import numpy as np
from fig5_chain import MODELS,model_inputs,read_cortex,make,artifact_endpoint,fit
from core import state
from protected_reml import _legacy_fit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
from panel_data import save,panel


def run(d,model,out,B,seed,design_labels):
    a,_=artifact_endpoint(make(d['X'],d['Y'],d),d['C'],d['pool'],'screen_r001')
    full,_,Y=fit(a,d['C'],a['k'],True)
    G=len(d['C']);ratios=np.arange(10,100,10)
    keep=np.array([j for j,label in enumerate(design_labels) if label=='intercept' or label.endswith(('_within','_between'))])
    if not len(keep):raise ValueError('No retained demographic design columns')
    Xgroup=d['X'][:,keep]
    group_full=_legacy_fit(state(Xgroup,Y,d['g'],d['C'],d['targets']),d['C'],0)[0]
    reference=np.stack([full,group_full])
    units=[np.flatnonzero(d['g']==j) for j in range(G)]
    correlations=np.full((2,len(ratios),B),np.nan)
    status=np.full(correlations.shape,'not_in_design',dtype='U20')
    for mode in range(2):
        for j,r in enumerate(ratios):
            if mode==1 and r<50:continue
            rng=np.random.RandomState(seed+(mode+1)*1000+r)
            for b in range(B):
                if mode==0:
                    rows=np.sort(np.concatenate([rng.choice(ix,max(1,len(ix)*r//100),replace=False) for ix in units]))
                else:
                    selected=rng.choice(G,max(2,G*r//100),replace=False)
                    rows=np.flatnonzero(np.isin(d['g'],selected))
                u,g=np.unique(d['g'][rows],return_inverse=True);c=d['C'][u];X=(d['X'] if mode==0 else Xgroup)[rows]
                design=np.c_[X,c[g]]
                if np.linalg.matrix_rank(design)<design.shape[1]:
                    status[mode,j,b]='rank_deficient';continue
                beta=_legacy_fit(state(X,Y[rows],g,c,d['targets']),c,0)[0]
                correlations[mode,j,b]=np.corrcoef(reference[mode],beta)[0,1]
                status[mode,j,b]='ok'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        summary=np.stack([np.nanmean(correlations,axis=2),np.nanquantile(correlations,.025,axis=2),np.nanquantile(correlations,.975,axis=2)],axis=-1)
    labels=['Subject','Unit'];fig='ExtendedDataFigure9'
    save(fig,out/model,dict(correlation_draws=correlations,status=status,seed=seed,proportions=ratios/100,series=summary,full_beta=full,reference_by_mode=reference,group_design_columns=keep,valid_draws=(status=='ok').sum(2)),[panel(fig,model+'_resampling','ribbons',x=ratios/100,series=summary,labels=np.array(labels,object))])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--harmonized',type=Path);p.add_argument('--models',choices=MODELS,nargs='+',default=list(MODELS))
    p.add_argument('--draws',type=int,default=100);p.add_argument('--seed',type=int,help='Override the model-specific fixed seed')
    args=p.parse_args()
    if args.draws<1:p.error('draws must be positive')
    for model in args.models:
        cohort=MODELS[model][0]
        if args.harmonized:
            with np.load(args.harmonized/(cohort+'_harmonized.npz')) as f:raw=dict(f)
        else:raw=read_cortex(args.inputs,cohort)
        d,_,_,labels=model_inputs(raw,args.inputs,model)
        seed=args.seed if args.seed is not None else {'PC1':2026101671,'PC2':2026101672,'USA_TL':2026101673}[model]
        run(d,model,args.output,args.draws,seed,labels)
