from pathlib import Path
import json
import numpy as np
import pandas as pd

T_WORDS = ['一致', '不准', '严格', '强制', '执行', '禁止', '约束', '维护', '遵从', '限制']

L_WORDS = ['允许', '创造', '变化', '开放', '灵活', '自主', '自由', '适应', '选择', '鼓励']

def weighted_score(x, polarity):
    (u, s, v) = np.linalg.svd(x - x.mean(0), full_matrices=False)
    orientation = np.sign(v[np.arange(len(v)), np.argmax(abs(v), axis=1)])
    scores = u * s * orientation
    loadings = v.T * orientation
    explained = s * s / (s * s).sum()
    k = int(np.searchsorted(np.cumsum(explained), 0.85) + 1)
    assert k == len(polarity)
    weights = explained[:k] / explained[:k].sum()
    score = scores[:, :k] * np.asarray(polarity) @ weights
    return (score, dict(retained_components=k, variance_threshold=0.85, explained_ratio=explained.tolist(), weights=weights.tolist(), component_polarity=polarity, loadings=loadings.tolist()))

def generate(word_file, province_file, output):
    words = T_WORDS + L_WORDS
    wide = pd.read_csv(word_file)[['names'] + words]
    province = pd.read_csv(province_file)
    merged = province[['names', 'z_tightness_control']].merge(wide, on='names', validate='one_to_one')
    assert len(merged) == 31
    matrix = merged[words].to_numpy(float)
    matrix /= matrix.sum(1, keepdims=True)
    (tight, tm) = weighted_score(matrix[:, :10], [1, 1, 1, 1])
    (loose, lm) = weighted_score(matrix[:, 10:], [1, -1])
    merged['pc1_T_clean'] = tight
    merged['pc1_L_clean'] = loose
    merged['pc1_Diff_clean'] = tight - loose
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output / 'model_b_cleaned_data.csv', index=False, encoding='utf-8-sig')
    (output / 'weighted_pca_parameters.json').write_text(json.dumps(dict(tightness=tm, looseness=lm), indent=2), encoding='utf-8')
    return merged
