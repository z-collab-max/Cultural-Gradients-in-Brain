import numpy as np

def fit_predict(x, y, z, model):
    if model == 'Additive':
        return z + np.mean(y - x)
    if model == 'Multiplicative':
        return z * np.dot(x, y) / np.dot(x, x)
    if model == 'ROI-specific additive':
        return z + y - x
    mu = x.mean()
    sd = x.std()
    u = (x - mu) / sd
    v = (z - mu) / sd

    def design(w):
        if model == 'Affine':
            return np.column_stack([np.ones(len(w)), w])
        if model.startswith('Polynomial'):
            return np.column_stack([w ** k for k in range(int(model[-1]) + 1)])
        centers = np.quantile(u, np.linspace(0, 1, 8))
        return np.column_stack([np.ones(len(w)), w, np.exp(-0.5 * ((w[:, None] - centers) / 0.5) ** 2)])
    a = design(u)
    b = design(v)
    penalty = np.eye(a.shape[1]) * len(x) * 0.001
    penalty[0, 0] = 0
    if model == 'Affine':
        penalty[:] = 0
    return b @ np.linalg.solve(a.T @ a + penalty, a.T @ y)
