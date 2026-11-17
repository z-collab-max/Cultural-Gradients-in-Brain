"""MATLAB REML fits using the shared China formula and canonical input fields."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import re
import numpy as np
from scipy.io import loadmat, savemat

INDIVIDUAL = ['age', 'sex', 'education'] + [f'genetic_pc{k}' for k in range(1, 11)] + [
    'father_education', 'mother_education', 'intracranial_volume', 'total_surface_area']
CONTROLS = INDIVIDUAL + [f'city_macro_pc{k}' for k in range(1,10)] + [f'province_macro_pc{k}' for k in range(1,3)]


def control_names(d):
    names=list(INDIVIDUAL)
    for prefix in ['city_macro_pc','province_macro_pc']:
        matched=sorted([k for k in d if re.fullmatch(prefix+r'[1-9][0-9]*',k)],key=lambda k:int(k[len(prefix):]))
        if not matched or matched!=[f'{prefix}{j+1}' for j in range(len(matched))]:
            raise ValueError('Supply consecutive PCs from separate macro preprocessing: '+prefix)
        names.extend(matched)
    return names

def z(x):
    x = np.asarray(x, dtype=float)
    sd = x.std(axis=0, ddof=1)
    if np.any(~np.isfinite(sd) | (sd <= 0)):
        raise ValueError('Every fitted column must have finite, positive variance')
    return (x - x.mean(axis=0)) / sd

def read_design(d):
    """Numeric group identifiers must encode birthplace province and city."""
    names=control_names(d)
    required = names + ['province_birth', 'province_city']
    missing = [k for k in required if k not in d]
    if missing:
        raise ValueError('Missing canonical China fields: ' + ', '.join(missing))
    x = np.column_stack([np.asarray(d[k]).ravel() for k in names])
    groups = np.column_stack([np.asarray(d[k]).ravel() for k in required[-2:]])
    if not np.isfinite(groups).all():
        raise ValueError('Province and city identifiers must be finite numeric codes')
    return x, groups

def fit_many(x, y, predictors, groups, names=None):
    """Return draw-by-ROI beta, t, p, df; refit variance components for every draw."""
    x, y, predictors, groups = map(np.asarray, (x, y, predictors, groups))
    if y.ndim == 1: y = y[:, None]
    if predictors.ndim == 1: predictors = predictors[:, None]
    names=list(CONTROLS if names is None else names)
    if names!=control_names(dict.fromkeys(names)) or x.shape != (len(y), len(names)) or groups.shape != (len(y), 2):
        raise ValueError('Expected named individual/macro-PC covariates and two group columns')
    if len(predictors) != len(y) or not all(np.isfinite(a).all() for a in (x,y,predictors,groups)):
        raise ValueError('Fit inputs must be finite and row-aligned')
    if any(len(np.unique(groups[:, k])) < 2 for k in range(2)):
        raise ValueError('Both random intercepts require at least two groups')
    z(x); z(y); z(predictors)
    executable = os.environ.get('MATLAB_EXECUTABLE') or shutil.which('matlab')
    if not executable:
        raise RuntimeError('MATLAB with Statistics and Machine Learning Toolbox is required; set MATLAB_EXECUTABLE')
    here = Path(__file__).resolve().parent
    quote = lambda p: str(p).replace("'", "''")
    with tempfile.TemporaryDirectory(prefix='china_mixed_') as tmp:
        source, target = Path(tmp)/'input.mat', Path(tmp)/'output.mat'
        savemat(source, dict(X=x,Y=y,predictors=predictors,groups=groups,covariate_names=np.array(names,dtype=object)))
        command = f"addpath('{quote(here)}'); study_china_batch('{quote(source)}','{quote(target)}')"
        subprocess.run([executable, '-batch', command], check=True, cwd=here)
        result = loadmat(target)
    return {k: result[k] for k in ('beta','t','p','df')}

def province_draws(values, groups, draws, seed, indices=None):
    if draws < 1:
        raise ValueError('draws must be positive')
    _, first, group = np.unique(groups[:, 0], return_index=True, return_inverse=True)
    scores = np.asarray(values)[first]
    if not np.allclose(scores[group], values):
        raise ValueError('Predictor must be constant within birthplace province')
    if indices is None:
        rng = np.random.default_rng(seed)
        indices = np.stack([rng.permutation(len(first)) for _ in range(draws)])
    indices = np.asarray(indices)
    if indices.shape != (draws, len(first)) or not np.issubdtype(indices.dtype, np.integer):
        raise ValueError('Expected draws-by-provinces integer permutation indices')
    if not np.all(np.sort(indices,axis=1) == np.arange(len(first))):
        raise ValueError('Each draw must permute every province exactly once')
    return np.column_stack([values, scores[indices].T[group]]), indices
