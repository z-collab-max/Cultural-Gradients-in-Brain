"""Equal-unit PCA with all nonzero components and reusable loadings."""
import numpy as np


def fit_pca(matrix, prefix):
    x = np.asarray(matrix, dtype=float)
    if x.ndim != 2 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError('PCA requires at least two complete geographic-unit rows')
    mean = x.mean(0)
    scale = x.std(0, ddof=1)
    scale[scale == 0] = 1
    z = (x - mean) / scale
    _, s, vt = np.linalg.svd(z, full_matrices=False)
    tol = max(z.shape) * np.finfo(float).eps * (s[0] if len(s) else 0)
    rank = int((s > tol).sum())
    loadings = vt[:rank].T.copy()
    for j in range(rank):
        if loadings[np.argmax(abs(loadings[:, j])), j] < 0:
            loadings[:, j] *= -1
    scores = z @ loadings
    names = [f'{prefix}_pc{j+1}' for j in range(rank)]
    total = float((s*s).sum())
    meta = dict(mean=mean.tolist(), scale=scale.tolist(), loadings=loadings.tolist(),
                singular_values=s.tolist(), rank=rank, tolerance=float(tol), names=names,
                explained_ratio=(s[:rank]**2/total).tolist() if total else [],
                units=len(x), weighting='equal geographic units')
    return scores, meta


def transform(matrix, metadata):
    x = np.asarray(matrix, float)
    if x.ndim != 2 or x.shape[1] != len(metadata['mean']) or not np.isfinite(x).all():
        raise ValueError('New macro inputs must match the fitted columns and be finite')
    return ((x-np.asarray(metadata['mean']))/np.asarray(metadata['scale'])) @ np.asarray(metadata['loadings'])
