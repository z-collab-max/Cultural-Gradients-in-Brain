import numpy as np

def _legacy_fit(s, c, k, return_y=False):
    cs = c[s['g']]
    x = s['X']
    y = s['Y']
    N = s['N']
    n = s['n']
    z = s['Z']
    D = np.c_[x, cs]
    dz = np.vstack([s['xz'], n * c])
    dd = np.block([[s['xx'], (x.T @ cs)[:, None]], [(cs @ x)[None, :], np.array([[cs @ cs]])]])
    dy = np.vstack([s['xy'], c @ s['zy']])
    zy = s['zy'].copy()
    yy = s['yy'].copy()
    if k:
        vs = cs / (s['se'] + n[s['g']] * s['sg'])
        E = cs - x @ np.linalg.solve(s['vxx'], x.T @ vs)
        (q, sv, _) = np.linalg.svd(np.c_[np.ones(N), cs, E], full_matrices=False)
        q = q[:, sv > sv[0] * 1e-12]
        gq = z.T @ q
        U = s['U'][:, :k, :]
        hmat = U @ gq
        rhs = np.einsum('rkg,gr->rk', U, s['zy']) - np.einsum('rkl,lr->rk', hmat, q.T @ y)
        gram = np.einsum('rkg,g,rlg->rkl', U, n, U) - hmat @ hmat.transpose(0, 2, 1)
        (eig, vec) = np.linalg.eigh((gram + gram.transpose(0, 2, 1)) / 2)
        inv = np.zeros_like(eig)
        np.divide(1, eig, out=inv, where=abs(eig) > 1e-10)
        coef = np.einsum('rki,ri->rk', vec, np.einsum('rki,rk->ri', vec, rhs) * inv)
        rm = np.einsum('rkg,rk->gr', U, coef)
        qc = np.einsum('rkl,rk->lr', hmat, coef)
        dy = dy - dz @ rm + D.T @ q @ qc
        zy = zy - n[:, None] * rm + gq @ qc
        yy = yy - np.sum(rhs * coef, axis=1)
        if return_y:
            yd = y - rm[s['g']] + q @ qc
    elif return_y:
        yd = y.copy()
    bols = np.linalg.solve(dd, dy)
    bz = np.linalg.solve(dd, dz)
    M = np.diag(n) - dz.T @ bz
    M = (M + M.T) / 2
    (eig, vec) = np.linalg.eigh(M)
    keep = eig > max(eig) * 1e-10
    eig = eig[keep]
    vec = vec[:, keep]
    a = vec.T @ (zy - dz.T @ bols)
    rss0 = yy - np.sum(dy * bols, axis=0)
    df = N - D.shape[1]

    def objective(lam):
        lam = np.asarray(lam)
        d = 1 + eig[:, None] * lam
        rss = rss0 - np.sum(a * a * lam / d, axis=0)
        assert np.all(rss > 0)
        return df * np.log(rss / df) + np.log(d).sum(0)
    grid = np.r_[0, np.logspace(-10, 6, 65)]
    values = np.array([objective(t) for t in grid])
    best = values.argmin(0)
    assert np.all(best < len(grid) - 1)
    lo = np.log(grid[np.maximum(1, best - 1)])
    hi = np.log(grid[np.minimum(len(grid) - 1, best + 1)])
    phi = (np.sqrt(5) - 1) / 2
    for _ in range(40):
        aa = hi - phi * (hi - lo)
        bb = lo + phi * (hi - lo)
        left = objective(np.exp(aa)) < objective(np.exp(bb))
        hi = np.where(left, bb, hi)
        lo = np.where(left, lo, aa)
    lam = np.exp((lo + hi) / 2)
    lam = np.where(objective(0) <= objective(lam), 0, lam)
    beta = bols[-1] - bz[-1] @ vec @ (a * lam / (1 + eig[:, None] * lam))
    sd = np.sqrt((yy - N * s['mean'] ** 2) / (N - 1))
    beta /= sd
    if return_y:
        return (beta, lam, yd)
    return (beta, lam)

def fit(s, c, k, return_y=False):
    if not k:
        return _legacy_fit(s, c, 0, return_y)
    from core import state
    X = s['X']
    g = s['g']
    n = s['n']
    cs = c[g]
    (_, anchor) = _legacy_fit(s, c, 0)
    Y = s['Y'].copy()
    for (j, l) in enumerate(anchor):
        w = l / (1 + n * l)
        vx = X - w[g, None] * s['xz'].T[g]
        vc = cs - w[g] * (n * c)[g]
        E = cs - X @ np.linalg.solve(X.T @ vx, X.T @ vc)
        h = E - w[g] * (s['Z'].T @ E)[g]
        (q, sv, _) = np.linalg.svd(np.c_[np.ones(len(cs)), h], full_matrices=False)
        q = q[:, sv > sv[0] * 1e-12]
        A = s['Z'] @ s['U'][j, :k].T
        A -= q @ (q.T @ A)
        gram = A.T @ A
        (ev, V) = np.linalg.eigh((gram + gram.T) / 2)
        iv = np.zeros_like(ev)
        np.divide(1, ev, out=iv, where=ev > 1e-10)
        Y[:, j] -= A @ (V @ (V.T @ (A.T @ s['Y'][:, j]) * iv))
    (b, lam) = _legacy_fit(state(X, Y, g, c, s['targets']), c, 0)
    return (b, lam, Y) if return_y else (b, lam)
