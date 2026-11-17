"""Convert saved simulation realizations into plotting arrays."""
from pathlib import Path
import argparse
import numpy as np
from scipy.io import loadmat, savemat
from scipy.stats import zscore


def export(source, output):
    output.mkdir(parents=True, exist_ok=True)
    for k in (1, 3, 5):
        s = loadmat(source/f'multi_artifact_K{k}.mat', simplify_cells=True)
        n = len(s['artif_rhos'])
        simulations = int(s['nSims_per_bin'])
        rois = s['maps_GT'].shape[1]
        assert s['maps_GT'].shape[0] == n*simulations
        for arm in ('M', 'C'):
            raw = s[f'corr_est_start_GT_{arm}']
            corrected = s[f'corr_est_end_GT_{arm}']
            ids = np.argmin(abs(raw-np.median(raw, axis=1)[:, None]), axis=1)
            indices = ids*n+np.arange(n)
            for short, reference in [('s', raw), ('e', corrected)]:
                a = zscore(s['maps_GT'], axis=1, ddof=1)
                b = zscore(s[f'maps_{short}{arm}'], axis=1, ddof=1)
                calculated = np.sum(a*b, axis=1)/(rois-1)
                np.testing.assert_allclose(calculated.reshape(simulations, n).T, reference, atol=1e-10)
            values = dict(artifact_rho=s['artif_rhos'], selected_simulation_1based=ids+1,
                          example_gt=zscore(s['maps_GT'][indices], axis=1, ddof=1),
                          example_raw=zscore(s[f'maps_s{arm}'][indices], axis=1, ddof=1),
                          example_corrected=zscore(s[f'maps_e{arm}'][indices], axis=1, ddof=1),
                          true_gradient_r=s['corr_GT_SA'], raw_gradient_r=s[f'corr_est_start_SA_{arm}'],
                          corrected_gradient_r=s[f'corr_est_end_SA_{arm}'], raw_fidelity=raw,
                          corrected_fidelity=corrected, source_count=k, arm=arm)
            savemat(output/f'SuppFig3_K{k}_{arm}.mat', values, do_compression=True)
            np.savez_compressed(output/f'SuppFig3_K{k}_{arm}.npz', **values)
            print(k, arm, 'PASS', flush=True)


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=root/'simulation/reference/si3/raw')
    parser.add_argument('--output', type=Path, default=root/'outputs/simulation_source/si3')
    args = parser.parse_args()
    export(args.source, args.output)
