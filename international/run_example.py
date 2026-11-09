"""Synthetic individual cortex -> within-unit ComBat -> de-artifacting -> REML."""
from pathlib import Path
import argparse
import json
import os
import sys

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parent/'core'))
import numpy as np
from scipy.io import savemat
from combat import combat_plan, combat
from model_spec import (z, individual_columns, NATIONS_DEMOGRAPHICS,
                        USA_DEMOGRAPHICS, NATIONS_BLOCKS, USA_MACROS)
from engine import make
from artifact import artifact_endpoint
from protected_reml import fit
from core import unit


def synthetic_cohort(groups, seed, n_rois=32):
    rng = np.random.default_rng(seed)
    g = np.repeat(np.arange(groups), 30)
    n = len(g)
    batch = 2*g + np.tile(np.repeat([0, 1], 15), groups)
    age = rng.uniform(18, 60, n)
    female = rng.integers(0, 2, n).astype(float)
    tiv = 1500 + 130*rng.normal(size=n)
    culture = z(rng.normal(size=(groups, 2)))
    axis = z(np.linspace(-1, 1, n_rois))
    raw = (2.5 + .08*culture[g, 0, None]*axis
           + .06*culture[g, 1, None]*np.sin(2*axis)
           + .002*(age[:, None]-35) + .03*female[:, None]
           + .0001*(tiv[:, None]-1500)
           + .08*rng.normal(size=(2*groups, n_rois))[batch]
           + .03*rng.normal(size=(groups, n_rois))[g]
           + .15*rng.normal(size=(n, n_rois)))
    return dict(g=g, batch=batch, age=age, female=female, tiv=tiv,
                culture=culture, raw=raw,
                macro=rng.normal(size=(groups, 12)),
                targets=np.column_stack([axis, z(rng.normal(size=n_rois))]),
                pool=np.array([rng.permutation(groups) for _ in range(10000)]))


def harmonize(d):
    a = dict(age=d['age'], sex_female_design=d['female'],
             sex_other_design=np.zeros(len(d['g'])),
             age_observed=np.ones(len(d['g']), bool),
             sex_observed=np.ones(len(d['g']), bool))
    plan = combat_plan(None, a, d['g'], d['batch'])
    return combat(d['raw'], plan)


def design(d, nations):
    demographics = individual_columns(d, nations)
    if nations:
        labels = list(NATIONS_DEMOGRAPHICS)
        columns = [demographics[key] for key in labels]
        offset = 0
        for block, names in NATIONS_BLOCKS.items():
            matrix = z(d['macro'][:, offset:offset+len(names)])
            u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
            scores = z(u*singular)
            columns.extend(scores[d['g'], j] for j in range(len(names)))
            labels.extend(f'{block}_PC{j+1}' for j in range(len(names)))
            offset += len(names)
    else:
        labels = list(USA_DEMOGRAPHICS)+list(USA_MACROS)
        columns = [demographics[key] for key in USA_DEMOGRAPHICS]
        columns.extend(z(d['macro'][d['g'], j]) for j in range(len(USA_MACROS)))
    X = np.column_stack(columns)
    assert np.isfinite(X).all()
    assert np.linalg.matrix_rank(np.c_[X, d['C'][d['g']]]) == X.shape[1]+1
    return X, labels


def permutation_test(state, d, observed, repetitions, seed):
    rng = np.random.default_rng(seed)
    permutations = np.array([rng.permutation(len(d['C'])) for _ in range(repetitions)])
    null_beta = np.empty((repetitions, state['Y'].shape[1]))
    null_r = np.empty((repetitions, d['targets'].shape[1]))
    counts = np.empty(repetitions, int)
    for j, order in enumerate(permutations):
        shuffled = d['C'][order]
        a, counts[j] = artifact_endpoint(state, shuffled, d['pool'], 'screen_r001')
        null_beta[j], _ = fit(a, shuffled, a['k'])
        null_r[j] = (unit(null_beta[j, :, None]).T@unit(d['targets']))[0]
    p = (abs(null_r) >= abs(observed)).mean(axis=0)
    return dict(permutations=permutations, null_beta=null_beta, null_r=null_r,
                p_two_sided_absolute=p, artifact_pool_counts=counts)


def run(output, permutations=20):
    output.mkdir(parents=True, exist_ok=True)
    summary = []
    for cohort, groups, seed in [('Nations', 26, 20261103), ('USA', 18, 20261104)]:
        raw = synthetic_cohort(groups, seed)
        np.savez_compressed(output/(cohort+'_synthetic_input.npz'), **raw)
        raw['harm'] = harmonize(raw)
        assert np.isfinite(raw['harm']).all()
        np.savez_compressed(output/(cohort+'_harmonized.npz'), harm=raw['harm'], g=raw['g'])
        models = [('PC1', 0), ('PC2', 1)] if cohort == 'Nations' else [('USA_TL', 0)]
        for model, column in models:
            d = dict(raw, C=raw['culture'][:, column])
            X, labels = design(d, cohort == 'Nations')
            Y = z(d['harm'])
            s = make(X, Y, d)
            for context, artifact in [('endpoint', artifact_endpoint)]:
                a, retained = artifact(s, d['C'], d['pool'], 'screen_r001', verify=True)
                beta, variance_ratio, cleaned = fit(a, d['C'], a['k'], True)
                correlation = (unit(beta[:, None]).T@unit(d['targets']))[0]
                assert np.isfinite(beta).all() and np.isfinite(cleaned).all()
                assert a['k'] == groups-1
                stem = model+'_'+context
                result = dict(beta=beta, variance_ratio=variance_ratio,
                              cleaned=cleaned, X=X, labels=np.array(labels),
                              culture=d['C'], g=d['g'], targets=d['targets'],
                              correlation=correlation, artifact_pc_count=a['k'],
                              artifact_cumulative_variance=a['cum_explained'])
                np.savez_compressed(output/(stem+'.npz'), **result)
                savemat(output/(stem+'.mat'), result)
                np.savetxt(output/(stem+'_source_data.csv'),
                           np.c_[np.arange(1, len(beta)+1), beta, d['targets']],
                           delimiter=',', header='roi,beta,synthetic_axis1,synthetic_axis2', comments='')
                row = dict(model=model, context=context, subjects=len(Y), groups=groups,
                           nuisance_columns=len(labels), labels=labels, k=a['k'],
                           screened_shuffles=retained, correlation=correlation.tolist())
                if permutations:
                    null = permutation_test(s, d, correlation, permutations, seed+column+100)
                    np.savez_compressed(output/(stem+'_permutation.npz'), **null)
                    savemat(output/(stem+'_permutation.mat'), null)
                    row['permutation_repetitions'] = permutations
                    row['permutation_p'] = null['p_two_sided_absolute'].tolist()
                summary.append(row)
                print(stem, 'PASS', 'k=', a['k'], 'X=', X.shape, flush=True)
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1]/'outputs/international')
    parser.add_argument('--permutations', type=int, default=20)
    args = parser.parse_args()
    if args.permutations < 0:
        parser.error('--permutations must be nonnegative')
    run(args.output.resolve(), args.permutations)
