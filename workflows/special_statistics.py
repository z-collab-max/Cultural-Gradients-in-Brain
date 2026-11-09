import numpy as np
from scipy.stats import rankdata


def cultural_rsa(cortex,culture,province,repetitions,seed,method):
    if method not in ['pearson','spearman']:raise ValueError('Specify pearson or spearman')
    groups=np.unique(province)
    brain=np.array([cortex[province==g].mean(0) for g in groups])
    scores=np.array([culture[province==g].mean(0) for g in groups])
    n=len(groups);i,j=np.triu_indices(n,1);brain_rdm=brain[i]-brain[j]
    def normalize(v):
        if method=='spearman':v=rankdata(v,axis=0)
        v=v-v.mean(0);return v/np.linalg.norm(v,axis=0)
    B=normalize(brain_rdm);observed=normalize(scores[i]-scores[j]).T@B
    rng=np.random.default_rng(seed);orders=np.array([rng.permutation(n) for _ in range(repetitions)])
    null=np.array([normalize(scores[o[i]]-scores[o[j]]).T@B for o in orders])
    p=(abs(null)>=abs(observed)).mean(0)
    return dict(effect_maps=observed,null_maps=null,rsa_p=p,province_mean_brain=brain,
                province_mean_culture=scores,brain_rdm=brain_rdm,culture_rdm=scores[i]-scores[j],orders=orders)


def cohen(a,b):
    a=np.asarray(a);b=np.asarray(b)
    pooled=((len(a)-1)*a.var(0,ddof=1)+(len(b)-1)*b.var(0,ddof=1))/(len(a)+len(b)-2)
    if np.any(pooled<=0):raise ValueError('Positive pooled variance required')
    return (a.mean(0)-b.mean(0))/np.sqrt(pooled)


def relocation_bootstrap(groups,repetitions,seed):
    # Each group is individuals x 2 columns: positive and negative ROI means.
    HH,LL,HL,LH,MHH,MLL=[np.asarray(groups[k]) for k in ['HH','LL','HL','LH','MHH','MLL']]
    for a in [HH,LL,HL,LH,MHH,MLL]:
        if a.ndim!=2 or a.shape[1]!=2 or len(a)<2 or not np.isfinite(a).all():raise ValueError('Each group requires at least two complete rows and two ROI-mean columns')
    n=min(len(HL),len(LH));rng=np.random.default_rng(seed);draws=[]
    for b in range(repetitions):
        hi=rng.integers(len(HH),size=n);lo=rng.integers(len(LL),size=n)
        draws.append(cohen(HH[hi],LL[lo]))
    draws=np.asarray(draws);difference=draws[:,0]-draws[:,1]
    cross=cohen(HL,LH);same=cohen(MHH,MLL);observed=np.array([cross[0]-cross[1],same[0]-same[1]])
    p=np.minimum(1,2*np.minimum((difference[:,None]<=observed).mean(0),(difference[:,None]>=observed).mean(0)))
    return dict(cohen_draws=draws,null_difference=difference,observed=observed,p=p,sample_size=n)


def traveling(raw1,raw2):
    from traveling_calibration import fit_predict
    raw=np.stack([raw1,raw2]);reference=np.median(raw,axis=1)
    if raw.ndim!=3 or raw.shape[1]<3 or not np.isfinite(raw).all() or np.any(raw<=0):raise ValueError('Aligned positive cortical matrices from at least three shared centers are required')
    deviation=100*(raw/reference[:,None,:]-1)
    spatial=np.array([np.corrcoef(deviation[0,j],deviation[1,j])[0,1] for j in range(raw.shape[1])])
    models=['Raw','ROI-specific additive','Additive','Multiplicative','Affine','Polynomial 2','Polynomial 3','RBF ridge']
    errors=np.empty((2,len(models),raw.shape[1],raw.shape[2]))
    for test in range(2):
        train=1-test
        for m,model in enumerate(models):
            for site in range(raw.shape[1]):
                ar=np.median(np.delete(raw[train],site,axis=0),axis=0)
                br=np.median(np.delete(raw[test],site,axis=0),axis=0)
                pred=raw[test,site] if model=='Raw' else fit_predict(raw[train,site],ar,raw[test,site],model)
                errors[test,m,site]=100*(pred-br)/br
    mae=np.median(abs(errors),axis=(2,3))
    return dict(heatmaps=deviation,spatial_r=spatial,median_r=np.median(spatial),
                relative_error=100*mae/mae[:,[0]],models=np.array(models),errors=errors,median_absolute_error=mae,
                span_medians=np.median(np.ptp(deviation,axis=2),axis=1))
