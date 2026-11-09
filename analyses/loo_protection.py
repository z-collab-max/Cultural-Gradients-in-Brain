import numpy as np

def directions(s, c, lam):
    X = s['X']
    g = s['g']
    n = s['n']
    cs = c[g]
    E = []
    H = []
    for l in lam:
        w = l / (1 + n * l)
        vx = X - w[g, None] * s['xz'].T[g]
        vc = cs - w[g] * (n * c)[g]
        e = cs - X @ np.linalg.solve(X.T @ vx, X.T @ vc)
        E.append(e)
        H.append(e - w[g] * (s['Z'].T @ e)[g])
    return (np.array(E).T, np.array(H).T)

def clean(s, H):
    Y = s['Y'].copy()
    ranks = []
    errors = []
    for j in range(Y.shape[1]):
        (q, sv, _) = np.linalg.svd(np.c_[np.ones(len(Y)), H[:, j]], full_matrices=False)
        q = q[:, sv > sv[0] * 1e-12]
        A = s['Z'] @ s['U'][j, :s['k']].T
        A -= q @ (q.T @ A)
        gram = A.T @ A
        (ev, V) = np.linalg.eigh((gram + gram.T) / 2)
        iv = np.zeros_like(ev)
        np.divide(1, ev, out=iv, where=ev > 1e-10)
        removed = A @ (V @ (V.T @ (A.T @ Y[:, j]) * iv))
        Y[:, j] -= removed
        ranks.append(int((ev > 1e-10).sum()))
        errors.append(float(np.max(abs(q.T @ removed))))
    return (Y, max(errors), sorted(set(ranks)))
