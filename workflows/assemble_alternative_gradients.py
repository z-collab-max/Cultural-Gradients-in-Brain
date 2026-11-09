"""Assemble six-gradient PC models, spatial spins and Webb intervals."""
import argparse
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analyses'))
from webb_statistics import cluster_score_intervals
from map_statistics import standardize
from panel_data import panel, save

NAMES=['SAaxis','Eigenmode2','Developmental','Evolutionary','Intersubject','CognitivePC1']

def tail(null,obs):
    return (abs(null)>=abs(obs)).mean(0)

def build(main,spatial,out,repetitions,seed):
    roi=np.asarray(main['roi_names']).astype(str).ravel()
    if not np.array_equal(roi,spatial['roi_names']):raise ValueError('ROI order differs')
    labels=list(spatial['gradient_labels']);G=standardize(spatial['gradients'][:,[labels.index(n) for n in NAMES]])
    spins=spatial['spin_indices']
    if spins.dtype.kind not in 'iu' or spins.ndim!=2 or spins.shape[1]!=len(roi) or spins.min()<0 or spins.max()>=len(roi):raise ValueError('Expected zero-based rotation by ROI indices')
    for pc,(sign,primary) in enumerate([(-1,1),(-1,0)]):
        point=main['beta'][pc]*sign;y=standardize(point)
        webb=cluster_score_intervals(point,main['cluster_influence'][pc]*sign,G,repetitions,seed,weights=main.get('webb_weights'))
        null=standardize(main['null_beta'][:,pc,:].T*sign).T
        perm=null@np.linalg.pinv(G).T;beta=webb['partial_beta']
        cols=np.arange(6) if pc==0 else np.array([0,2,3,4,5])
        spin_maps=standardize(point[spins].T).T
        projector=np.linalg.pinv(G[:,cols]);sb=projector@y;sn=spin_maps@projector.T
        r=webb['r'][primary];rn=null@G[:,primary]/(len(y)-1);rs=spin_maps@G[:,primary]/(len(y)-1)
        source=dict(roi_names=roi,gradient_names=np.array(NAMES),spin_gradient_names=np.asarray(NAMES)[cols],
                    standardized_gradients=G,point=point,partial_beta=beta,province_null=perm,p_province=tail(perm,beta),
                    spin_beta=sb,spin_null=sn,p_spin=tail(sn,sb),r=r,zero_order_province_null=rn,zero_order_spin_null=rs,
                    zero_order_p_province=tail(rn,r),zero_order_p_spin=tail(rs,r),webb_r_ci=webb['r_ci'],webb_beta_ci=webb['beta_ci'])
        panels=[]
        for test,observed,draws,names in [('province',beta,perm,NAMES),('spin',sb,sn,np.asarray(NAMES)[cols])]:
            panels.append(panel('ExtendedDataFigure14',f'PC{pc+1}_{test}','violin',beta=observed,null=draws,p=tail(draws,observed),labels=np.asarray(names,object)))
        pidx=int(np.flatnonzero(cols==primary)[0])
        panels.append(panel('ExtendedDataFigure14',f'PC{pc+1}_gradient','scatter',x=G[:,primary],y=y,lh=np.char.startswith(roi,'lh_'),
                            xlabel=NAMES[primary],ylabel='Culture effect',stat_null=sn[:,pidx],stat_observed=sb[pidx],stats=f'β = {beta[primary]:.3f}; p = {tail(sn,sb)[pidx]:.3f}'))
        ci=webb['beta_ci']
        panels.append(panel('ExtendedDataFigure14',f'PC{pc+1}_Webb','intervals',estimate=beta,low=ci[0],high=ci[1],labels=np.array(NAMES,object),xlabel='Coefficient (beta)'))
        save('ExtendedDataFigure14',out/f'PC{pc+1}',source,panels)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['main','spatial','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--draws',type=int,default=10000);p.add_argument('--seed',type=int,default=42)
    a=p.parse_args()
    with np.load(a.spatial) as f:build(loadmat(a.main,simplify_cells=True),dict(f),a.output,a.draws,a.seed)

