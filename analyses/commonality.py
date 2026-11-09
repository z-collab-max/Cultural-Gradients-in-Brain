import numpy as np

def regress(y, columns):
    X = np.column_stack([np.ones(len(y))] + columns)
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ coef
    return (1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum(), coef, pred)
