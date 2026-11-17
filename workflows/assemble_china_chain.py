"""Build Fig. 2 gradient and Fig. 3 commonality sources from a shared China run."""
import argparse
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from assemble import assemble
from panel_data import save
from run_specialized import run


def build(results, spatial_file, output, repetitions, seed):
    main = loadmat(results / 'main_maps.mat', simplify_cells=True)
    with np.load(spatial_file, allow_pickle=False) as archive:
        spatial = dict(archive)
    roi_names = np.asarray(main['roi_names']).astype(str).ravel()
    if not np.array_equal(roi_names, spatial['roi_names'].astype(str).ravel()):
        raise ValueError('Spatial data and model maps require the same ordered ROI names')
    data = dict(spatial, effect_maps=main['beta'], null_maps=main['null_beta'],
                cluster_influence=main['cluster_influence'])
    source, panels = assemble('Figure2', data, repetitions, seed)
    source.update(roi_names=roi_names, p_model=main['p_model'], p_permutation_roi=main['p_permutation'])
    save('Figure2', output, source, panels)
    maps = loadmat(results / 'commonality.mat', simplify_cells=True)
    run('commonality', {k:v for k,v in maps.items() if not k.startswith('__')}, output, repetitions, seed)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--spatial', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    build(args.results, args.spatial, args.output, args.repetitions, args.seed)
