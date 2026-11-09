"""Generate China's two cultural-PC trajectories from individual residuals."""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse,sys
from pathlib import Path
import json
import h5py
import numpy as np
from scipy.io import savemat, loadmat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from panel_data import panel as plot_panel,save as save_panel

from china_mixed import CONTROLS, z, read_design, fit_many

def unit(a):
    a = a - a.mean(axis=-1, keepdims=True)
    return a / np.linalg.norm(a, axis=-1, keepdims=True)

def generate(brain_file, raw_file, culture_file, index_file, target_file, output_dir, draws=2000):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with h5py.File(raw_file) as h:
        valid = h['valididx'][()].ravel().astype(bool)
        levels, pcs = h['raw_res'].shape
    assert (levels, pcs) == (29, 2)
    with h5py.File(brain_file) as h:
        d = h['brain_data']
        controls,groups=read_design(d);controls=controls[valid];groups=groups[valid]
        real = np.column_stack([d[f'pc7_{i}'][()].ravel()[valid] for i in (1, 2)])
    if h5py.is_hdf5(culture_file):
        with h5py.File(culture_file) as h:
            culture = np.column_stack([h[f'Culdata/pc7_{i}'][()].ravel() for i in (1, 2)])
    else:
        d = loadmat(culture_file, simplify_cells=True)['Culdata']
        culture = np.column_stack([np.asarray(d[f'pc7_{i}']).ravel() for i in (1, 2)])
    with h5py.File(target_file) as h:
        targets = unit(np.stack([-h['mode2'][()].ravel(), -h['SAaxis'][()].ravel()]))
    indices = np.load(index_file)['indices_zero_based'][:, :draws]
    _, first, group = np.unique(groups[:, 0], return_index=True, return_inverse=True)
    assert indices.shape == (len(first), draws)
    assert np.array_equal(real, real[first][group])
    n = len(real)
    G = np.eye(len(first))[group]
    counts = G.sum(axis=0)
    metadata = dict(brain_file=str(brain_file), raw_file=str(raw_file),
                    culture_file=str(culture_file), index_file=str(index_file),
                    target_file=str(target_file), controls=CONTROLS,
                    participants=n, permutations=draws,
                    raw_res_columns=list(range(1, 30)),
                    removed_PCs=[0] + list(range(2, 30)),
                    permuted_predictor='direct province shuffle followed by sample-SD standardization')
    for pc in range(2):
        if not np.allclose(culture,real[first]):raise ValueError('Culture rows must match sorted province identifiers')
        if not np.all(np.sort(indices,axis=0)==np.arange(len(first))[:,None]):raise ValueError('Invalid province permutations')
        predictors=np.column_stack([real[:,pc],culture[indices,pc][group]])
        null_r = np.empty((draws, levels))
        observed_r = np.empty(levels)
        observed_beta = np.empty((levels, len(targets[pc])))
        null_beta = np.lib.format.open_memmap(out / f'PC{pc+1}_null_beta.npy', mode='w+',
                          dtype='float64', shape=(levels, draws, len(targets[pc])))
        with h5py.File(raw_file) as h:
            for k in range(levels):
                y = h[h['raw_res'][k, pc]][()].T
                result=fit_many(controls,y,predictors,groups)
                beta=result['beta'][1:];observed_beta[k]=result['beta'][0]
                null_beta[k] = beta
                null_r[:, k] = unit(beta) @ targets[pc]
                observed_r[k] = unit(observed_beta[k]) @ targets[pc]
        null_beta.flush()
        p = np.mean(np.abs(null_r) >= np.abs(observed_r), axis=0)
        name = ['China_PC1_vs_eigenmode2', 'China_PC2_vs_SA'][pc]
        panel = dict(name=name, r=observed_r, p=p, selected=28,
                     components=np.array([0]+list(range(2,30))),
                     **{f'null_{k}': null_r[:, k] for k in range(levels)})
        stem = f'EDFig2d_{pc+1}'
        np.savez_compressed(out / (stem + '.npz'), **panel)
        savemat(out / (stem + '.mat'), panel, do_compression=True)
        save_panel('ExtendedDataFigure2',out/f'PC{pc+1}',panel,[plot_panel('ExtendedDataFigure2',name,'trajectory',**{k:v for k,v in panel.items() if k!='name'})])
        np.save(out / f'PC{pc+1}_observed_beta.npy', observed_beta)
        print(f'PC{pc+1}: completed {levels} levels x {draws} permutations', flush=True)
    (out / 'generation_inputs.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('brain_file', 'raw_file', 'culture_file', 'index_file', 'target_file', 'output_dir'):
        p.add_argument('--' + key.replace('_', '-'), required=True)
    p.add_argument('--draws',type=int,default=2000)
    generate(**vars(p.parse_args()))
