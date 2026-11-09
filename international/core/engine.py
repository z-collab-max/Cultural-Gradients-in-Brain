import numpy as np
from core import prepare_metric

def make(X, Y, d):
    s = prepare_metric(X, Y, d['g'], d['C'], d['targets'])
    Q = np.linalg.qr(X, mode='reduced')[0]
    RZ = s['Z'] - Q @ (Q.T @ s['Z'])
    RY = Y - Q @ (Q.T @ Y)
    K = RZ.T @ RZ
    T = RZ.T @ RY
    (ev, V) = np.linalg.eigh((K + K.T) / 2)
    ok = ev > max(ev) * 1e-10
    ev = ev[ok]
    V = V[:, ok]
    F = V * np.sqrt(ev)
    tt = V.T @ T / np.sqrt(ev)[:, None]
    (ii, jj) = np.triu_indices(len(ev))
    s.update(K=K, T=T, F=F, tri_i=ii, tri_j=jj, tp=tt[ii] * tt[jj])
    return s
