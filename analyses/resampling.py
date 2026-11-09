import numpy as np


def wild_cluster_draws(fitted,residuals,clusters,statistic,repetitions=10000,seed=20261103,weights='webb'):
    fitted=np.asarray(fitted);residuals=np.asarray(residuals)
    if fitted.shape!=residuals.shape:raise ValueError('Fitted and residual shapes differ')
    _,group=np.unique(clusters,return_inverse=True);rng=np.random.default_rng(seed)
    if weights=='webb':support=np.array([-np.sqrt(1.5),-1,-np.sqrt(.5),np.sqrt(.5),1,np.sqrt(1.5)])
    elif weights=='rademacher':support=np.array([-1,1])
    else:raise ValueError('Unknown weights')
    draws=[]
    for b in range(repetitions):
        w=rng.choice(support,group.max()+1)[group]
        if residuals.ndim==2:w=w[:,None]
        draws.append(statistic(fitted+residuals*w))
    return np.asarray(draws)


def snr_curve(signal,residuals,groups,fit_map,levels,repetitions=100,seed=20261103):
    """SNR levels are signal/noise mean-square power ratios (variance ratios for centered inputs)."""
    signal=np.asarray(signal,float);residuals=np.asarray(residuals,float)
    if signal.shape!=residuals.shape:raise ValueError('Signal and residual shapes differ')
    rng=np.random.default_rng(seed);groups=np.asarray(groups);rows=[np.flatnonzero(groups==g) for g in np.unique(groups)]
    noise_scale=np.sqrt(np.mean(residuals**2));signal_scale=np.sqrt(np.mean(signal**2))
    if signal_scale==0 or noise_scale==0:raise ValueError('Nonzero signal and noise required')
    maps=[]
    for level in levels:
        if level<=0:raise ValueError('SNR levels must be positive')
        current=[]
        for b in range(repetitions):
            noise=np.empty_like(residuals)
            for ix in rows:noise[ix]=residuals[rng.choice(ix,len(ix),replace=True)]
            y=signal+noise*(signal_scale/(np.sqrt(level)*noise_scale))
            current.append(fit_map(y))
        maps.append(current)
    return np.asarray(maps)
