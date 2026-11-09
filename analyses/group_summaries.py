import numpy as np


def summaries(values,groups):
    values=np.asarray(values);groups=np.asarray(groups)
    if values.ndim==1:values=values[:,None]
    labels=np.unique(groups);mean=[];sem=[];n=[]
    for g in labels:
        a=values[groups==g];n.append(len(a));mean.append(a.mean(0));sem.append(a.std(0,ddof=1)/np.sqrt(len(a)))
    return dict(labels=labels,mean=np.array(mean),sem=np.array(sem),n=np.array(n))


def shared_roi_bootstrap(positive,negative,groups,repetitions=10000,seed=20261103):
    positive=np.asarray(positive);negative=np.asarray(negative);groups=np.asarray(groups)
    labels=np.unique(groups);rows=[np.flatnonzero(groups==g) for g in labels]
    count=min(map(len,rows));rng=np.random.default_rng(seed)
    draws=np.empty((repetitions,len(labels),2))
    for b in range(repetitions):
        for j,indices in enumerate(rows):
            selected=rng.choice(indices,count,replace=True)
            draws[b,j]=[positive[selected].mean(),negative[selected].mean()]
    return dict(groups=labels,draws=draws,sample_size=count,
                interval=np.quantile(draws,[.025,.975],axis=0))


def grid_counts(xy,cell_size):
    # xy contains projected coordinates in the same distance unit as cell_size.
    cells=np.floor(np.asarray(xy)/cell_size).astype(np.int64)
    labels,counts=np.unique(cells,axis=0,return_counts=True)
    return labels,counts


def grid_overlap(xy_a,xy_b,cell_size):
    a,na=grid_counts(xy_a,cell_size);b,nb=grid_counts(xy_b,cell_size)
    aa={tuple(k):int(n) for k,n in zip(a,na)};bb={tuple(k):int(n) for k,n in zip(b,nb)}
    cells=sorted(set(aa)|set(bb))
    return np.array(cells),np.array([[aa.get(k,0),bb.get(k,0)] for k in cells])
