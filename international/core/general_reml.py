import numpy as np

def fit(s, c, k, return_y=False):
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
        E = vs - s['vx'] @ np.linalg.solve(s['vxx'], x.T @ vs)
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
