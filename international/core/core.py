import numpy as np
from general_reml import fit

def z(a):
    return (a - a.mean(0)) / a.std(0, ddof=1)

def unit(a):
    a = a - a.mean(0)
    return a / np.linalg.norm(a, axis=0)

def state(X, Y, g, C, targets, sg=1.0, se=1.0):
    n = np.bincount(g)
    Z = np.eye(len(n))[g]
    xz = X.T @ Z
    den = se + n * sg
    vx = (X - (sg / den)[g, None] * xz.T[g]) / se
    return dict(X=X, Y=Y, g=g, gid=g, C=C, N=len(g), n=n, ng=n, Z=Z, xz=xz, xx=X.T @ X, xy=X.T @ Y, zy=Z.T @ Y, yy=np.sum(Y * Y, 0), mean=Y.mean(0), vx=vx, vxx=X.T @ vx, sg=sg, se=se, targets=targets)

def prepare_metric(X, Y, g, C, targets):
    a = state(X[:, 1:], Y, g, C, targets)
    cc = np.ones(len(a['n']))
    (_, lam) = fit(a, cc, 0)
    s = state(X, Y, g, C, targets)
    bols = np.linalg.solve(s['xx'], s['xy'])
    bz = np.linalg.solve(s['xx'], s['xz'])
    M = np.diag(s['n']) - s['xz'].T @ bz
    (ev, V) = np.linalg.eigh((M + M.T) / 2)
    ok = ev > max(ev) * 1e-10
    ev = ev[ok]
    V = V[:, ok]
    t = V.T @ (s['zy'] - s['xz'].T @ bols)
    rss = s['yy'] - np.sum(s['xy'] * bols, 0) - np.sum(t * t * lam / (1 + ev[:, None] * lam), 0)
    se = rss / (len(g) - X.shape[1])
    s = state(X, Y, g, C, targets, float(np.median(lam * se)), float(np.median(se)))
    B = s['xz'].T / (s['se'] + s['n'] * s['sg'])[:, None]
    A = s['vxx']
    K = np.diag(s['n'] / (s['se'] + s['n'] * s['sg'])) - B @ np.linalg.solve(A, B.T)
    T = s['zy'] / (s['se'] + s['n'] * s['sg'])[:, None] - B @ np.linalg.solve(A, s['vx'].T @ Y)
    (ke, kv) = np.linalg.eigh((K + K.T) / 2)
    good = ke > max(ke) * 1e-10
    F = kv[:, good] * np.sqrt(ke[good])
    tt = kv[:, good].T @ T / np.sqrt(ke[good])[:, None]
    assert np.max(abs(F @ tt - T)) < 1e-07
    (i, j) = np.triu_indices(sum(good))
    s.update(K=K, T=T, F=F, tri_i=i, tri_j=j, tp=tt[i] * tt[j])
    return s
