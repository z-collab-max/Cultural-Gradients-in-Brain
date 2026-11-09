import numpy as np
from core import unit

def artifact_endpoint(s, c, pool, mode, verify=False):
    assert mode == 'screen_r001', 'Expected screen_r001 artifact scores'
    S = c[pool].T
    corr = unit(S).T @ unit(c[:, None])
    S = S[:, abs(corr[:, 0]) < 0.01]
    count = S.shape[1]
    if count < 2:
        raise ValueError('fewer than2 pool shuffles satisfy abs(r)<0.01')
    n = s['n']
    G = len(n)
    S -= (n @ S / s['N'])[None, :]
    S /= np.sqrt(n @ (S * S) / (s['N'] - 1))[None, :]
    den = np.sum(S * (s['K'] @ S), axis=0)
    if min(den) <= 1e-10:
        raise ValueError('degenerate screened score')
    beta = S.T @ s['T'] / den[:, None]
    W = np.sqrt(n[:, None]) * S
    grams = np.einsum('gp,rp,hp->rgh', W, beta.T ** 2, W, optimize=True)
    if verify:
        for r in sorted(set([0, s['Y'].shape[1] // 2, s['Y'].shape[1] - 1])):
            A = W * beta[:, r]
            assert np.max(abs(A @ A.T - grams[r])) < 1e-06
    (ev, V) = np.linalg.eigh((grams + grams.transpose(0, 2, 1)) / 2)
    ev = np.maximum(ev[:, ::-1], 0)
    V = V[:, :, ::-1]
    expl = 100 * ev[:, :G - 1] / ev.sum(1)[:, None]
    k = G - 1
    if not 1 <= k <= G - 1:
        raise ValueError('invalid artifact PC count')
    a = s.copy()
    a.update(U=V[:, :, :G - 1].transpose(0, 2, 1) / np.sqrt(n)[None, None, :], k=k, cum_explained=np.cumsum(expl.mean(0)))
    return (a, count)
