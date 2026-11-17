from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

T_WORDS = ['一致', '不准', '严格', '强制', '执行', '禁止', '约束', '维护', '遵从', '限制']

L_WORDS = ['允许', '创造', '变化', '开放', '灵活', '自主', '自由', '适应', '选择', '鼓励']

def rank_normal_scores(scores):
    """Column-wise Blom transform using average ranks for ties."""
    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 2 or scores.shape[0] < 2 or not np.isfinite(scores).all():
        raise ValueError('Rank normalization requires a finite observations-by-scores matrix')
    ranks = rankdata(scores, method='average', axis=0)
    return norm.ppf((ranks - 0.375) / (len(scores) + 0.25))

def pc1_score(x):
    """Centered PCA; orient each loading by its largest absolute coefficient."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or min(x.shape) < 2 or not np.isfinite(x).all():
        raise ValueError('PCA requires a finite matrix with at least two rows and columns')
    (u, s, v) = np.linalg.svd(x - x.mean(0), full_matrices=False)
    if s[0] <= 0:
        raise ValueError('PCA requires nonconstant keyword frequencies')
    orientation = np.sign(v[np.arange(len(v)), np.argmax(abs(v), axis=1)])
    scores = u * s * orientation
    loadings = v.T * orientation
    explained = s * s / (s * s).sum()
    return scores[:, 0], dict(retained_components=1, explained_ratio=explained.tolist(),
        center=x.mean(0).tolist(), loadings=loadings.tolist(),
        orientation='largest_absolute_loading_positive')

def generate(word_file, province_file, output):
    words = T_WORDS + L_WORDS
    wide = pd.read_csv(word_file)[['names'] + words]
    province = pd.read_csv(province_file)
    merged = province[['names', 'z_tightness_control']].merge(wide, on='names', validate='one_to_one')
    assert len(merged) == 31
    matrix = merged[words].to_numpy(float)
    if not np.isfinite(matrix).all() or np.any(matrix < 0) or np.any(matrix.sum(1) <= 0):
        raise ValueError('Keyword counts must be finite, nonnegative and have positive row totals')
    matrix /= matrix.sum(1, keepdims=True)
    (tight, tm) = pc1_score(matrix[:, :10])
    (loose, lm) = pc1_score(matrix[:, 10:])
    merged['pc1_T_clean'] = tight
    merged['pc1_L_clean'] = loose
    merged['pc1_Diff_clean'] = tight - loose
    transformed = rank_normal_scores(np.column_stack([tight-loose, tight, loose]))
    merged[['rank_normal_Diff', 'rank_normal_T', 'rank_normal_L']] = transformed
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output / 'model_b_cleaned_data.csv', index=False, encoding='utf-8-sig')
    (output / 'pca_parameters.json').write_text(json.dumps(dict(tightness=tm, looseness=lm,
        rank_normalization=dict(method='Blom',offset=0.375,ties='average',
        score_order=['PC1_difference','tightness_PC1','looseness_PC1'])), indent=2), encoding='utf-8')
    return merged
