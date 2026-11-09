import argparse
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import rice_fields as model
from china_mixed import read_design, fit_many, province_draws
from scipy.io import savemat

parser = argparse.ArgumentParser()
parser.add_argument('--inputs', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
parser.add_argument('--draws', type=int, default=10000)
args = parser.parse_args()
cfg = json.loads(args.inputs.read_text(encoding='utf-8-sig'))
for key in cfg:
    cfg[key] = args.inputs.parent / cfg[key]
OUT = args.output
OUT.mkdir(parents=True, exist_ok=True)
N_PERM = args.draws
if N_PERM < 1:
    raise ValueError('draws must be positive')
with h5py.File(cfg['roi_order']) as f:
    names = model.matlab_cellstr(f, f['Talhelm_Gradient_Results/thickness/brain_var'])
subjects = pd.read_csv(cfg['subjects'], dtype=str).set_index('subject_id', verify_integrity=True)
indicators = pd.read_csv(cfg['indicators']).set_index('Province_CN', verify_integrity=True)
with h5py.File(cfg['brain']) as f:
    d = f['brain_data']
    ids = model.matlab_cellstr(f, d['MRIID']) if h5py.check_dtype(ref=d['MRIID'].dtype) else [str(int(x)) for x in d['MRIID'][()].ravel()]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate MRI IDs')
    provinces_all = subjects.loc[ids, 'province'].fillna('').str.strip().to_numpy()
    X,groups=read_design(d)
    raw = {k:model.numeric(d,k) for k in ['centerID_rerange']}
    brain = np.column_stack([model.numeric(d,k) for k in names])
matched = indicators.reindex(provinces_all)
rice = {field:matched[col].to_numpy(float) for field,col in zip(model.RICE_MODELS.values(),['PercentPaddy','EnvironmentalRiceSuitability'])}
valid = np.isfinite(np.column_stack([brain,X,*raw.values(),*rice.values()])).all(axis=1)
valid &= provinces_all != ''
prov = provinces_all[valid]
province_names, index = np.unique(prov, return_inverse=True)
n = len(index)
env = dict(names=names, brain=brain, rice=rice)
with h5py.File(cfg['roi_order'], "r") as f:
    area_names = model.matlab_cellstr(f, f["Talhelm_Gradient_Results/area/brain_var"])
with h5py.File(cfg['brain'], "r") as f:
    area = np.column_stack([model.numeric(f["brain_data"], k) for k in area_names])
assert np.isfinite(area[valid]).all()

assert [k.removesuffix("_area") for k in area_names] == [k.removesuffix("_thickness") for k in env["names"]]

with h5py.File(cfg['gradients'], "r") as f:
    gradients = np.column_stack([model.numeric(f, k)[:len(names)] for k in model.MAP_FIELDS])
gz = model.zcols(gradients)
with h5py.File(cfg['spins'], "r") as f:
    spins = np.asarray(f["perm_id"][()], int)[:N_PERM] - 1
assert spins.shape == (N_PERM, len(names)) and spins.min() >= 0 and spins.max() < len(names)
spun_gradients = [model.zrows(gradients[:, j][spins]) for j in range(6)]
projector = np.linalg.pinv(gz)

centers, cidx = np.unique(raw["centerID_rerange"][valid], return_inverse=True)
center_dummy = np.eye(len(centers))[cidx, :-1]
center_design = np.column_stack([np.ones(n), center_dummy])
assert np.linalg.matrix_rank(center_design) == center_design.shape[1]
rows, diagnostics, maps, plot_arrays = [], [], {}, {}


def p_empirical(null, observed):
    return float(np.mean(np.abs(null) >= abs(observed)))


def expcs_residualize(y_raw):
    """Match glmfit(K-1 centre dummies, ROI) followed by zscore(residual)."""
    residual = y_raw - center_design @ np.linalg.lstsq(center_design, y_raw, rcond=None)[0]
    residual = model.zcols(residual)
    # Residuals must be orthogonal to the intercept and every retained centre dummy.
    assert np.max(np.abs(center_design.T @ residual)) < 1e-8
    return residual


for phenotype, y_all in [("thickness", env["brain"]), ("area", area)]:
    y = expcs_residualize(y_all[valid])
    for label, field in model.RICE_MODELS.items():
        x=env['rice'][field][valid]
        candidates,_=province_draws(x,groups[valid],N_PERM,20260905)
        result=fit_many(X[valid],y,candidates,groups[valid]);beta=result['beta'][0]
        plot_arrays[f"{phenotype}_{label}_tVal"] = result['t'][0]
        plot_arrays[f"{phenotype}_{label}_pVal"] = result['p'][0]
        plot_arrays[f"{phenotype}_{label}_dfe"] = result['df'][0]
        perm_z=model.zrows(result['beta'][1:])
        bz = model.zvec(beta)
        direct_real = bz @ gz / (len(gz)-1)
        direct_province_null = perm_z @ gz / (len(gz)-1)
        real_multi = projector @ bz
        multi_province_null = perm_z @ projector.T
        multi_spin_null = model.zrows(beta[spins]) @ projector.T
        plot_arrays[f"{phenotype}_{label}_province_null"] = multi_province_null
        plot_arrays[f"{phenotype}_{label}_spin_null"] = multi_spin_null
        plot_arrays[f"{phenotype}_{label}_observed"] = real_multi
        for j, gradient in enumerate(model.MAP_NAMES):
            ps = p_empirical(spun_gradients[j] @ bz / (len(gz)-1), direct_real[j])
            pp = p_empirical(direct_province_null[:, j], direct_real[j])
            mps = p_empirical(multi_spin_null[:, j], real_multi[j])
            mpp = p_empirical(multi_province_null[:, j], real_multi[j])
            rows.append(dict(
                phenotype=phenotype, rice_indicator=label, gradient=gradient,
                r=float(direct_real[j]), p_spin=ps, p_province=pp,
                pass_both=bool(ps < .05 and pp < .05),
                multivariable_beta=float(real_multi[j]),
                multivariable_p_spin=mps, multivariable_p_province=mpp,
                multivariable_pass_both=bool(mps < .05 and mpp < .05)))
        assert np.isfinite(beta).all() and np.isfinite(perm_z).all()
        maps[f"{phenotype}_{label}"] = beta.tolist()
        plot_arrays[f"{phenotype}_{label}_map"] = beta

for phenotype in ('thickness','area'):
    for indicator in model.RICE_MODELS:
        prefix = phenotype+'_'+indicator
        v = {k[len(prefix)+1:]:x for k,x in plot_arrays.items() if k.startswith(prefix+'_')}
        rr = [r for r in rows if r['phenotype']==phenotype and r['rice_indicator']==indicator]
        for test in ('spin','province'):
            v['p_'+test] = np.array([r['multivariable_p_'+test] for r in rr])
        p = v['pVal']
        ps = np.sort(p)
        passes = ps <= .05*np.arange(1,len(ps)+1)/len(ps)
        v['fdr_significant'] = p<=ps[passes].max() if passes.any() else np.zeros(len(p),bool)
        v.update(roi_names=np.asarray(names),gradient_names=np.asarray(model.MAP_NAMES))
        np.savez_compressed(OUT/f'EDFig8_{prefix}.npz',**v)
        savemat(OUT/f'EDFig8_{prefix}.mat',{k:x.astype(object) if isinstance(x,np.ndarray) and x.dtype.kind=='U' else x for k,x in v.items()},do_compression=True)
print(f'{n} participants, {len(province_names)} provinces, {N_PERM} draws')
