import numpy as np

def combat_plan(meta, a, g, sid):
    P = np.c_[a['age'].ravel(), a['sex_female_design'].ravel(), a['sex_other_design'].ravel()]
    seen = np.c_[a['age_observed'].ravel(), a['sex_observed'].ravel(), a['sex_observed'].ravel()].astype(bool)
    plan = []
    for j in np.unique(g):
        ix = np.flatnonzero(g == j)
        (_, b) = np.unique(sid[ix], return_inverse=True)
        nb = b.max() + 1
        B = np.eye(nb)[b]
        protect = np.empty((len(ix), 0))
        active = []
        if nb > 1:
            for k in range(3):
                rows = seen[ix][:, active + [k]].all(1)
                D = np.c_[B[rows], protect[rows]]
                if np.linalg.matrix_rank(D) == D.shape[1] and np.linalg.matrix_rank(np.c_[D, P[ix][rows, k]]) == D.shape[1] + 1:
                    protect = np.c_[protect, P[ix, k]]
                    active.append(k)
        plan.append(dict(ix=ix, b=b, B=B, protect=protect, active=active))
    return plan

def combat(Y, plan):
    out = Y.copy()
    for p in plan:
        ix = p['ix']
        B = p['B']
        nb = B.shape[1]
        if nb == 1:
            continue
        y = Y[ix]
        D = np.c_[B, p['protect']]
        n = B.sum(0)
        N = len(ix)
        coef = np.linalg.solve(D.T @ D, D.T @ y)
        grand = n / N @ coef[:nb]
        vp = np.mean((y - D @ coef) ** 2, axis=0)
        vp[vp == 0] = np.median(vp[vp != 0])
        stand = np.tile(grand, (N, 1)) + D[:, nb:] @ coef[nb:]
        sd = (y - stand) / np.sqrt(vp)
        gh = B.T @ sd / n[:, None]
        dh = np.array([sd[p['b'] == j].var(0, ddof=1) for j in range(nb)])
        gb = gh.mean(1)
        t2 = gh.var(1, ddof=1)
        dm = dh.mean(1)
        dv = dh.var(1, ddof=1)
        ap = (2 * dv + dm ** 2) / dv
        bp = (dm * dv + dm ** 3) / dv
        adj = sd.copy()
        for j in range(nb):
            hit = p['b'] == j
            go = gh[j].copy()
            do = dh[j].copy()
            for iteration in range(10000):
                gn = (t2[j] * n[j] * gh[j] + do * gb[j]) / (t2[j] * n[j] + do)
                dn = (0.5 * np.sum((sd[hit] - gn) ** 2, axis=0) + bp[j]) / (n[j] / 2 + ap[j] - 1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    change = max(np.max(abs(gn - go) / go), np.max(abs(dn - do) / do))
                (go, do) = (gn, dn)
                if change <= 0.001:
                    break
            else:
                raise RuntimeError('ComBat did not converge')
            adj[hit] = (sd[hit] - gn) / np.sqrt(dn)
        out[ix] = adj * np.sqrt(vp) + stand
    assert np.isfinite(out).all()
    return out
