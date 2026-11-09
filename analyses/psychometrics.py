import numpy as np
import pandas as pd
from scipy import stats

def anova_icc(y, groups):
    d = pd.DataFrame({'y': pd.to_numeric(y, errors='coerce'), 'group': groups}).dropna()
    if d['group'].dtype == object:
        d = d[d.group.astype(str).str.len() > 0]
    samples = [v.y.to_numpy() for (_, v) in d.groupby('group')]
    n = len(d)
    k = len(samples)
    mu = d.y.mean()
    msb = sum((len(v) * (v.mean() - mu) ** 2 for v in samples)) / (k - 1)
    msw = sum((((v - v.mean()) ** 2).sum() for v in samples)) / (n - k)
    (F, p) = stats.f_oneway(*samples)
    icc1 = (msb - msw) / (msb + (n / k - 1) * msw)
    return {'n': n, 'groups': k, 'F': float(F), 'p': float(p), 'ICC1': float(icc1), 'ICC2': float((msb - msw) / msb)}

def rwg(data, items, maximum):
    null = (maximum ** 2 - 1) / 12
    rows = []
    for (name, g) in data.dropna(subset=['from_province']).groupby('from_province'):
        if len(g) < 2:
            continue
        mv = g[items].apply(pd.to_numeric, errors='coerce').var(ddof=1).mean()
        j = len(items)
        v = 0.0 if mv > null else j * (1 - mv / null) / (j * (1 - mv / null) + mv / null)
        rows.append({'province': str(name), 'n': len(g), 'item_variance': float(mv), 'rwg': float(v)})
    return (float(np.mean([v['rwg'] for v in rows])), rows)
