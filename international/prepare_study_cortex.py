"""Assemble individual cortical measurements and harmonize within cultural units."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def demographics(config):
    m = pd.read_csv(config['nations_demographics'])
    u = pd.read_csv(config['usa_demographics'])
    assert m.original_row.is_unique and u.subid.is_unique
    stage = read_json(config['demographic_changes'])
    decisions = read_json(config['demographic_decisions'])
    changes = list(stage['verified_changes'])
    exclusions = list(decisions['exclusions'])
    for x in read_json(config['pnas_identity_resolution'])['subjects']:
        sid = 'OpenNeuro_'+x['dataset']+'_'+x['participant']
        ix = u.index[u.subid.eq(sid)]
        assert len(ix) == 1
        i = int(ix[0])
        if not x['eligible_18_60']:
            exclusions.append(dict(cohort='USA', index=i))
        for field, old, new in [
            ('age', x['old_age'], x['source_age']),
            ('sex', 1 if x['old_sex'] == 'M' else 2,
             1 if x['source_sex'] == 'M' else 2)]:
            if old != new:
                changes.append(dict(cohort='USA', index=i, field=field, old=old, new=new))
    for x in changes:
        frame = m if x['cohort'] == 'Nations' else u
        assert abs(float(frame.at[x['index'], x['field']])-x['old']) < 1e-7
        frame.at[x['index'], x['field']] = x['new']
    for x in decisions['age_midpoint_flags']:
        m.at[x['index'], 'age_estimated'] = True
    result = {}
    for cohort, frame in [('Nations', m), ('USA', u)]:
        keep = np.ones(len(frame), bool)
        keep[[x['index'] for x in exclusions if x['cohort'] == cohort]] = False
        result[cohort] = frame.loc[keep].copy()
    return result


def build(config_path, output):
    config_path = Path(config_path).resolve()
    config = {k: str(config_path.parent / v) for k, v in read_json(config_path).items()}
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    demo = demographics(config)
    with h5py.File(config['nations_cortical_table']) as f:
        group = f['nationsdata_all']
        field_meta = f[group.attrs['MATLAB_fields']][()]
        fields = [b''.join(x).decode() for x in field_meta]
        roi = [x.removesuffix('_thickness') for x in fields
               if x.endswith('_thickness') and '500_aparc' not in x and 'unknown' not in x]
        assert len(roi) == len(set(roi)) == 308
        ids = [''.join(chr(c) for c in f[x][()].ravel())
               for x in group['lh_500_aparc_thickness'][()].ravel()]
        values = np.column_stack([group[x+'_thickness'][()].ravel() for x in roi])
    keep = [i for i, sid in enumerate(ids) if not sid.startswith('Israel-')
            and not (sid.startswith('Germany-3-') and 'acq-MP2R' in sid)]
    nids = [ids[i].replace('Turkey-1-', 'Mexico-3-').replace('Spain-2-', 'Belgium-3-') for i in keep]
    blocks = [values[keep]]
    usa_append = loadmat(config['nations_appended_usa_table'], simplify_cells=True)
    columns = [list(usa_append['roi_names']).index(x) for x in roi]
    for j, state in enumerate(['Louisiana', 'Missouri', 'Wisconsin'], 1):
        ix = np.flatnonzero(usa_append['state'] == state)
        blocks.append(usa_append['thickness'][ix][:, columns])
        nids += [f'USA-{j}-'+str(usa_append['subid'][i]) for i in ix]
    nids = np.array(nids)
    nm = demo['Nations']
    ni = nm.original_row.to_numpy(int)-1
    assert np.array_equal(nids[ni], nm.participant.to_numpy())
    nraw = np.vstack(blocks)[ni]
    usa = pd.concat([pd.read_csv(config['usa_cortical_table']),
                     pd.read_csv(config['usa_batch3_cortical_table'])], ignore_index=True)
    assert usa.subid.is_unique
    usa = usa.set_index('subid')
    um = demo['USA']
    selected_usa = usa.loc[um.subid]
    assert np.array_equal(selected_usa.state, um.state)
    assert np.array_equal(selected_usa.dataset, um.dataset)
    uraw = selected_usa[roi].to_numpy(float)
    sys.path.insert(0, str(Path(__file__).resolve().parent/'core'))
    from combat import combat, combat_plan
    summary = []
    for cohort, frame, raw, group_col, batch_col in [
        ('Nations', nm, nraw, 'country', 'corrected_batch'),
        ('USA', um, uraw, 'state', 'dataset')]:
        is_nations = cohort == 'Nations'
        g, group_names = pd.factorize(frame[group_col], sort=False)
        sid, batch_names = pd.factorize(frame[batch_col], sort=False)
        age = frame['age_updated' if is_nations else 'age'].to_numpy(float)
        female = frame.female.to_numpy(float) if is_nations else frame.sex.eq(2).to_numpy(float)
        age_observed = (~frame.age_estimated).to_numpy(float) if is_nations else np.ones(len(frame))
        sex_observed = (~frame.sex_estimated).to_numpy(float) if is_nations else np.ones(len(frame))
        data = dict(raw=raw, g=g, sid=sid, age=age, sex_female_design=female,
                    sex_other_design=np.zeros(len(frame)), age_observed=age_observed,
                    sex_observed=sex_observed, female=female, tiv=frame.eTIV.to_numpy(float),
                    roi_names=np.array(roi), participant=frame['participant' if is_nations else 'subid'].to_numpy(str),
                    group_names=np.array(group_names, str), batch_names=np.array(batch_names, str),
                    original_input_index=frame.index.to_numpy(int))
        assert np.isfinite(raw).all() and np.isfinite(data['tiv']).all()
        assert age.min() >= 18 and age.max() <= 60
        plan = combat_plan(None, data, g, sid)
        data['harm'] = combat(raw, plan)
        np.savez_compressed(output/(cohort+'_harmonized.npz'), **data)
        audit = [dict(group=str(group_names[g[p['ix'][0]]]), n=len(p['ix']),
                      batches=p['B'].shape[1], protected_columns=p['active']) for p in plan]
        (output/(cohort+'_combat_plan.json')).write_text(json.dumps(audit, indent=2), encoding='utf-8')
        summary.append(dict(cohort=cohort, n=len(frame), rois=len(roi), groups=len(group_names)))
        print(cohort, raw.shape, 'ComBat saved', flush=True)
    paths = {k: v for k, v in config.items() if isinstance(v, str) and Path(v).is_file()}
    provenance = dict(inputs={k: dict(sha256=hashlib.sha256(Path(v).read_bytes()).hexdigest())
                             for k, v in paths.items()}, cohorts=summary)
    (output/'input_hashes.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.config.resolve(), args.output.resolve())
