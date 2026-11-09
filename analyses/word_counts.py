from pathlib import Path
import json
import numpy as np
import pandas as pd
from word_scores import T_WORDS,L_WORDS

def generate(raw_file, province_file, output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    raw = pd.read_excel(raw_file, sheet_name='Sheet1')
    d = raw.copy()
    province = '省/直辖市'
    word = '关键词'
    date = '日期'
    value = '出现关键词的微博指数'
    words = T_WORDS + L_WORDS
    shanxi = d[province].str.contains('山西', na=False)
    changes = []
    for (term, rule) in [('限制', 'date'), ('自由', 'threshold')]:
        select = shanxi & d[word].eq(term)
        extreme = d[date].eq(20190531) if rule == 'date' else d[value].gt(100000)
        mask = select & extreme
        normal = select & ~extreme
        replacement = float(d.loc[normal, value].median())
        if not mask.any() or not np.isfinite(replacement):
            raise ValueError(f'Cannot apply {term} replacement')
        for index in d.index[mask]:
            changes.append({'province': str(d.loc[index, province]), 'keyword': term, 'date': str(d.loc[index, date]), 'old_index': float(d.loc[index, value]), 'replacement': replacement, 'normal_dates': int(normal.sum())})
        d.loc[mask, value] = replacement
    before = raw.groupby([province, word])[value].sum().unstack().fillna(0)[words]
    after = d.groupby([province, word])[value].sum().unstack().fillna(0)[words]
    observed = pd.read_csv(province_file)[['names', 'z_tightness_control']]

    def match(s):
        matches = [name for name in observed.names if s in name]
        if len(matches) > 1:
            raise ValueError(f'Province mapping is not unique: {s}')
        return matches[0] if matches else None
    after['names'] = [match(s) for s in after.index]
    after = after.dropna(subset=['names'])
    before = before.loc[after.index]
    wide = observed.merge(after.reset_index(drop=True), on='names', validate='one_to_one')
    assert len(wide) == 31
    wide.to_csv(out / 'cleaned_word_counts.csv', index=False, encoding='utf-8-sig')
    raw_wide = before.copy()
    raw_wide['names'] = [match(s) for s in before.index]
    raw_wide.reset_index(drop=True).to_csv(out / 'raw_word_counts.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(changes).to_csv(out / 'daily_replacements.csv', index=False, encoding='utf-8-sig')
    comparison = before.to_numpy() != after[words].to_numpy()
    (out / 'cleaning_parameters.json').write_text(json.dumps({'raw_file': str(Path(raw_file).resolve()), 'province_file': str(Path(province_file).resolve()), 'specified_date': 20190531, 'freedom_threshold': 100000, 'replacements': changes, 'aggregate_cells_changed': int(comparison.sum()), 'aggregate_cells': int(comparison.size)}, ensure_ascii=False, indent=2), encoding='utf8')
    return out / 'cleaned_word_counts.csv'
