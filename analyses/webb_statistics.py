"""Cluster-score Webb draws and spatial-statistic percentile intervals."""
import numpy as np


def z(x,axis=0):
    x=np.asarray(x,float);sd=x.std(axis=axis,ddof=1,keepdims=True)
    if not np.isfinite(x).all() or np.any(sd<=0):raise ValueError('Finite nonconstant maps required')
    return (x-x.mean(axis=axis,keepdims=True))/sd


def balanced_webb(repetitions,groups,seed=42):
    if repetitions<2 or repetitions%2:raise ValueError('An even repetition count >=2 is required')
    support=np.array([-np.sqrt(1.5),-1,-np.sqrt(.5),np.sqrt(.5),1,np.sqrt(1.5)])
    rng=np.random.RandomState(seed);base=np.resize(support,repetitions//2)
    half=np.column_stack([base[rng.permutation(len(base))] for _ in range(groups)])
    return np.concatenate([half,-half])


def interval_summary(draws,confidence=.95,quantile_method='hazen'):
    a=np.asarray(draws,float)
    if len(a)<2 or not np.isfinite(a).all():raise ValueError('At least two finite draws required')
    tail=(1-confidence)/2
    return dict(ci=np.quantile(a,[tail,1-tail],axis=0,method=quantile_method),
                se=a.std(axis=0,ddof=1),mean=a.mean(axis=0),
                sign_tail_p=np.minimum(1,2*np.minimum((a<=0).mean(0),(a>=0).mean(0))))


def spatial_statistics(point,draws,gradients,orientation=1,quantile_method='hazen'):
    if orientation not in [-1,1]:raise ValueError('Orientation must be +1 or -1')
    point=np.asarray(point).ravel()*orientation;draws=np.asarray(draws)*orientation
    G=z(gradients);p=z(point);B=z(draws,axis=1)
    if G.shape[0]!=len(p) or B.shape[1]!=len(p):raise ValueError('ROI dimensions differ')
    if np.linalg.matrix_rank(G)!=G.shape[1]:raise ValueError('Gradient design is rank deficient')
    projector=np.linalg.solve(G.T@G,G.T)
    r=p@G/(len(p)-1);rd=B@G/(len(p)-1)
    beta=projector@p;bd=B@projector.T
    out=dict(point_map=point,r=r,partial_beta=beta,r_draws=rd,beta_draws=bd)
    for prefix,values in [('r',rd),('beta',bd)]:
        for key,value in interval_summary(values,quantile_method=quantile_method).items():out[prefix+'_'+key]=value
    return out


def cluster_score_intervals(point,influence,gradients,repetitions=10000,seed=42,orientation=1,weights=None):
    point=np.asarray(point).ravel();influence=np.asarray(influence,float)
    if influence.ndim!=2 or influence.shape[1]!=len(point):raise ValueError('Influence must be groups by ROIs')
    if weights is None:weights=balanced_webb(repetitions,len(influence),seed)
    weights=np.asarray(weights,float)
    if weights.ndim!=2 or weights.shape[1]!=len(influence):raise ValueError('Weights must be draws by groups')
    maps=point[None,:]+weights@influence
    out=spatial_statistics(point,maps,gradients,orientation)
    out.update(weights=weights,map_draws=maps,seed=seed)
    return out
