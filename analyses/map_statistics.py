import numpy as np
from scipy.stats import rankdata


def standardize(x):
    x = np.asarray(x, float)
    scale = x.std(axis=0, ddof=1)
    if np.any(scale == 0) or not np.isfinite(x).all():
        raise ValueError('Expected finite, nonconstant columns')
    return (x-x.mean(axis=0))/scale


def map_correlation(x, y, rank=False):
    x, y = np.asarray(x).ravel(), np.asarray(y).ravel()
    if rank:
        x, y = rankdata(x), rankdata(y)
    return float(standardize(x)@standardize(y)/(len(x)-1))


def spin_test(x, y, spin_indices):
    # spin_indices: rotations x ROIs, zero-based atlas correspondences.
    x, y = np.asarray(x).ravel(), np.asarray(y).ravel()
    indices = np.asarray(spin_indices)
    if indices.ndim != 2 or indices.shape[1] != len(x):
        raise ValueError('Spin array must be rotations by ROIs')
    if not np.issubdtype(indices.dtype, np.integer) or indices.min()<0 or indices.max()>=len(x):
        raise ValueError('Invalid zero-based spin indices')
    observed = map_correlation(x,y)
    null = np.array([map_correlation(x[i],y) for i in indices])
    return observed, null, np.mean(abs(null)>=abs(observed))


def partial_map_coefficients(effect, gradients):
    y=standardize(np.asarray(effect).ravel());X=standardize(gradients)
    design=np.column_stack([np.ones(len(y)),X])
    if np.linalg.matrix_rank(design)!=design.shape[1]:
        raise ValueError('Gradient design is rank deficient')
    return np.linalg.lstsq(design,y,rcond=None)[0][1:]


def commonality(target, deculturation, acculturation, relocation):
    from commonality import regress
    total,coef,predicted=regress(target,[deculturation,acculturation,relocation])
    ci=regress(target,[deculturation])[0];ca=regress(target,[acculturation])[0]
    both=regress(target,[deculturation,acculturation])[0]
    parts=np.array([both-ca,both-ci,ci+ca-both,total-both,total])
    assert np.isclose(parts[:4].sum(),total)
    return parts,coef,predicted


def bh_fdr(p):
    p=np.asarray(p,float);order=np.argsort(p);q=np.empty_like(p)
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q
