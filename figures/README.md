# Figure rendering

Install requirements.txt from the repository root. MATLAB statistical layouts require Statistics and Machine Learning Toolbox. Geography exports SVG/PNG; all MATLAB layouts save editable FIG files. The shared style is Arial 8 pt and 0.85 pt line width.

```sh
python figures/run_layout_demo.py --output outputs/layout_demo
```

```matlab
addpath('figures'); addpath('figures/layouts');
p=fullfile(pwd,'outputs','layout_demo'); out=fullfile(p,'rendered');
render_brain_map(fullfile(p,'brain.mat'),fullfile(out,'brain'));
render_vertex_map(fullfile(p,'vertex.mat'),fullfile(out,'vertex'));
reproduce_suppfig1(fullfile(p,'si1.mat'),out);
reproduce_suppfig2(fullfile(p,'si2.mat'),out);
reproduce_suppfig4(fullfile(p,'si4.mat'),out);
plot_fig2e_snr(fullfile(p,'snr.mat'),out);
reproduce_edfig15a(fullfile(p,'ed15a.mat'),out);
reproduce_edfig15b(fullfile(p,'ed15b.mat'),out);
```

## Cortical maps

render_brain_map accepts a MAT file with roi_names (cell strings such as lh_precentral_part1_thickness), values [R], color_limits [2], and optional significant [R] logical mask. It matches each name to the bundled atlas and renders lateral/medial views of both hemispheres. Omitted atlas regions are gray. Thus a 307-region vector is supported when its ROI-name list identifies the omitted region. Significance contours use the original midpoint-triangle boundary algorithm.

For ED4, render_vertex_map accepts values_lh [10242], values_rh [10242] and color_limits [2] in the bundled fsaverage5 vertex order. Corresponding meshes and cortical masks are included. Vertex preprocessing and model fitting remain upstream operations.

## Geographic maps

```sh
python figures/render_maps.py --input inputs/maps.json --output outputs/maps
```

The JSON contains columns and panels. Each panel has region (china, world or usa), title, optional extent [west,east,south,north], and points. Each point contains longitude, latitude, count and optionally composition (nonnegative counts for pie sectors) or color. Inputs are city/site aggregates. Optional region_values maps exact boundary names to scalar culture values, with color_limits [low,high]. World/USA polygon fills support these values; China files supply boundary lines and city markers. Names can be inspected in assets/geography/*.json. The six-panel demo illustrates center/origin maps. Region-specific inset placement and final manuscript panel lettering remain caller-configurable.

## Dedicated statistical layouts

All nine functions in layouts accept (source_file, output_dir). input_schemas.json lists required MAT field names/types and reference dimensions; dimensions are illustrative where sample size is variable. No values are embedded. The layouts are:

| Entry | Content |
|---|---|
| plot_fig2a | Cultural PCA biplot and explained variance; Fig2a_Data input; display flips PC2 |
| plot_fig2e_snr | Two-PC SNR recovery curves, theory and calibrated uncertainty |
| plot_fig2e_replacement | Prioritized/random replacement curves, lower-tail probabilities and Holm correction; source fields from `workflows/assemble_fig2_details.py` |
| plot_fig3jk_ed11b | Positive/negative ROI bootstrap distributions and ROI interval panels; Fig4_ED9_Data and colorpair inputs |
| reproduce_edfig15a | Cultural-space PCA with region labels, loadings and group ellipses |
| reproduce_edfig15b | Joint/local/world PCA axes; display flips joint PC1 |
| reproduce_suppfig1 | Word dictionary, variance and external-validity fits |
| reproduce_suppfig2 | Traveling-participant deviation maps and calibration comparisons |
| reproduce_suppfig4 | Relocation duration and the most frequent cross-province routes |

For SI4, workflows/assemble.py computes the required fields from duration and origin/destination labels; pass the resulting source_data.mat to the dedicated plotter. For SI2, the same workflow produces layout fields for aligned shared centers; retaining extra nonshared centers requires separate descriptive heatmap rows. The joint cultural PCA generator is analyses/cultural_pca.py:joint_cultural_pca; supply gender-balanced province means, aligned country dimensions, the gender-balanced overall China sample mean and the China world-reference score. Sign orientation is fixed before plotting. Regression and bootstrap results supplied to plotters must come from the intended calculation branch; plotting does not recompute them.
