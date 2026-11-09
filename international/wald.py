import numpy as np


def wald_statistics(d, beta, lam, clean):
    g = d['g']
    design = np.c_[d['X'], d['C'][g]]
    n = np.bincount(g)
    sx = np.array([design[g == i].sum(0) for i in range(len(n))])
    sy = np.array([clean[g == i].sum(0) for i in range(len(n))])
    se, t = [], []
    for j, ratio in enumerate(lam):
        weight = ratio/(1+n*ratio)
        a = design.T@design-sx.T@(weight[:, None]*sx)
        v = design.T@clean[:, j]-sx.T@(weight*sy[:, j])
        coef = np.linalg.solve(a, v)
        residual = clean[:, j]-design@coef
        grouped = np.bincount(g, weights=residual)
        variance = (residual@residual-np.sum(weight*grouped**2))/(len(g)-design.shape[1])
        stderr = np.sqrt(variance*np.linalg.inv(a)[-1, -1])
        scale = clean[:, j].std(ddof=1)
        assert abs(coef[-1]/scale-beta[j]) < 1e-6
        se.append(stderr/scale)
        t.append(coef[-1]/stderr)
    return np.array(se), np.array(t)


