# Specialized source-data workflows

Run from the repository root after installing `requirements.txt`. These examples generate synthetic inputs. Supply your own aligned measurements and fitted results for empirical analyses.

```sh
python workflows/run_remaining_demo.py --output outputs/specialized
python workflows/run_specialized.py --analysis gradients --input inputs/gradients.mat --output outputs/source
```

Each analysis writes a figure directory with numeric source files and `render_manifest.json`. Execute all plotting jobs in that directory with:

```matlab
addpath('figures');
render_workflow('outputs/source/gradients/ExtendedDataFigure3','outputs/figures/gradients');
```

Statistical plots save editable `.fig` files and PNG previews. Optional `roi_names` must match the bundled atlas; this additionally exports surface maps for gradients, hierarchy RSA/slope maps, effect PCA and fitted-map branches. All ROI arrays and spatial rotations must have identical ordering. Demo permutations of ROI indices test interfaces only; empirical spin tests require valid spatial rotations.

## Input fields

N denotes observations, G groups, R ROIs, M models, B draws, T tasks and D cultural dimensions. Indices in NPZ inputs are zero-based. Named labels must be Unicode arrays. MAT inputs are also accepted; use NPZ to preserve singleton dimensions. Except for optional fields explicitly named below, all listed fields are required. `run_remaining_demo.py` produces one NPZ per route and an `input_schema.json` with example dimensions.

| Analysis | Figure | Input fields and shape |
|---|---|---|
| province_summary | Fig. 2 | `x`, `y` [N,M], `province` [N]; means and SEM across individuals within province |
| ecological_comparison | ED1 | Same as province_summary; optional `x_label`, `y_label` |
| demographics | ED1 | `continuous` [N,V], `continuous_labels` [V], `group` [N]; optional `categorical` [N] |
| psychometrics | ED1 | `dimension_scores` [N,D], `items` [N,I], `item_mask` [D,I] Boolean, `anova_group`, `icc_group` [N], `response_max`, `compute_alpha`, `dimension_labels` [D]; alpha is omitted where compute_alpha is false |
| brain_icc | ED1 | `adjusted_brain` [M,N,R], `province` [N], `measure_labels` [M]; model adjustment must precede assembly |
| gradients | ED3 | Nine [R] vectors: `SAaxis`, `mode2_display`, `xu_evo`, `reardon_dev`, `intervar`, `cogPC1`, `myelin`, `gene_PC1`, `fcgradient1`; alternatively use original `DK308_modes_asym` and `fcgradientsM`, taking negative mode 2 and FC gradient 1 |
| commonality | Fig. 3 | `mover_maps`, `stayer_maps`, `deculturation_maps`, `acculturation_maps`, `relocation_maps` [M,R]; target is mover minus stayer |
| hierarchy | Fig. 4 / ED12 | `task_effects` [branches,R,T], `branch_labels`, `task_labels`, `example_roi`, `spin_indices` [B,R], `left_hemisphere` [R]; selects one best-fitting order separately for each branch |
| rice_culture | ED8 | `group_values` [G,V], `labels` [V] containing rice and culture variables |
| effect_pca | ED15 | `dimension_effects` [R,D], `culture_effects` [R,M], `reference_component` [M], `dimension_labels` [D], `spin_indices` [B,R], `left_hemisphere` [R] |
| resampling | ED9 | `full_maps` [M,R], `refit_maps` [M,B,P,R], `proportions` [P], `series_labels` [M] |
| network | ED11 | `HH`, `LL`, `HL`, `LH`, `MHH`, `MLL` [group N,2], paired positive/negative ROI means; this network panel uses the smaller strict bootstrap tail with a 1/B floor; Fig. 3 bootstrap uses the separate two-sided route in assemble.py |
| reference_histograms | ED6 | `reference` [R], `effect_maps` [M,R], `null_maps` [B,M,R] |
| dimension_omission | Additional sensitivity | `full_maps` [M,R], `fold_maps` [D,M,R], `dimension_labels` [D], `spin_indices` [B,R]; match and orient components using cultural scores upstream |
| thresholds | ED2 | `endpoint_map` [R], `threshold_maps` [P,R], `threshold_labels` [P] |
| overlap_curves | ED9 | `scales` [P], `overlap_counts` [P,3], ordered both / imaging-only / culture-only |
| combat_versions | ED5 | Fitted-map contract below, with explicit configuration labels |
| area_volume | ED7 | Fitted-map contract; labels distinguish measure and predictor |
| rice_brain | ED8 | Fitted-map contract; labels distinguish rice share/suitability and cortical measure |
| international | Fig. 5 | Fitted-map contract plus `webb_map_draws` [B,M,R] from full endpoint refits |
| word_layout | SI1 | `word_counts` [G,20], `scores` [G,3] ordered difference/tightness/looseness, `validity_score` [G], `word_labels` [20], `dictionary` [10,4], `province_names` [G] |
| traveling_unshared | SI2 | `participant1` [C1,R], `participant2` [C2,R], `centers1` [C1], `centers2` [C2]; positive thickness; calibration uses common centers, descriptive heatmaps also retain unshared centers |

Fitted-map contract: `effect_maps` [M,R], `null_maps` [B,M,R], `gradients` [R,K], `spin_indices` [S,R], `model_labels` [M]; optional `gradient_labels` [K], `left_hemisphere` [R], `roi_names` [R]. Outer null maps must come from unrestricted shuffles. Webb draws and permutation nulls are distinct inputs.

## Model and preprocessing connections

`china/run_model_branches.m` accepts a MAT struct array `branches` with `name`, `T`, `cortex` [N,R], `predictor_names` and optional `additional_covariates`. Each branch uses the common China covariate and random-intercept formula. `T` follows `china/INPUT_SCHEMA.md`; extra predictors must be numeric columns. It saves `beta`, `se`, `p` and `model_labels` separately for every branch.

```matlab
addpath('china');
run_model_branches('inputs/branches.mat','outputs/models');
run_dimension_omission('inputs/dimension_omission.mat','outputs/dimension_omission.mat');
generate_gradient_inputs('inputs/gradient_summaries','outputs/gradients.mat');
```

Dimension-omission inputs are `T`, `cortex`, `cultural_dimensions` [G,D], `province_ids` [G] and optional `dimension_labels`. Recomputed cultural PCs are matched and oriented against full cultural PCA scores before cortical fitting. Add valid `spin_indices` to its saved result for `dimension_omission` assembly.

The original gradient builder reads `DK308_summary_L.mat` / `DK308_summary_R.mat` containing `Summary_DK308_L` / `Summary_DK308_R`, and `surface_500_L_164k_midthickness_summary.mat` / `surface_500_R_164k_midthickness_summary.mat` containing `surface_500_L` / `surface_500_R`. Hemisphere summaries require SAaxis, intervar, cogPC1, xu_evo, reardon_dev, xu_FChomology, gene_PC1, T1w_T2w_T1wT2w and fcgradient. Surface structures require label, modes and cortex. Label zero is excluded.

`prepare_combat_versions.py` takes `cortex` [N,R], `batch` [N], `biological_covariates` [N,C] without intercept, `culture_pc1`, `culture_pc2`, `recruit` [N]. It saves the PC1/PC2-protected and PC1/PC2-plus-recruitment-protected matrices separately. Fit each matrix using its intended downstream model, then pass the resulting maps to assembly.

```sh
python workflows/prepare_combat_versions.py --input inputs/combat.npz --output outputs/combat
python workflows/bridge_china_results.py --results outputs/china/source_data/hierarchy_inputs.mat --spatial inputs/hierarchy_spatial.npz --output outputs/source
python workflows/bridge_word_scores.py --input outputs/words/model_b_cleaned_data.csv --output outputs/source
```

The China bridge reads `mediation` and `behavior_order` from `hierarchy_inputs.mat`, written in MATLAB v7 format by `export_china_source` and `run_china_example`. Spatial input contains `spin_indices`, `left_hemisphere`, `example_roi` and optional `roi_names`. It reconstructs original task labels, then selects one best order for behavior and each cultural mediation branch. `analyses/run_word_pipeline.py` also calls the word bridge automatically after cleaning and scoring. `international/dimension_loo.py` provides endpoint cultural-dimension refits using the same individual design structure as `international/run_example.py`.

## Collect separately saved map results

```sh
python workflows/collect_map_results.py --manifest inputs/models.json --output outputs/source
```

Paths in this manifest are relative to the manifest file:

```json
{
  "analysis": "international",
  "spatial": "spatial.npz",
  "models": [
    {"label": "Nations_PC1", "effect_file": "PC1_endpoint.npz", "null_file": "PC1_endpoint_permutation.npz", "webb_file": "PC1_webb.npz"}
  ]
}
```

Default keys are `beta` [R], `null_beta` [B,R] and `map_draws` [B,R]. Override them with `effect_key`, `null_key`, `webb_key`. For combined arrays `beta` [M,R] and nulls [B,M,R], specify zero-based `model_index`. Shared spatial inputs supply gradients, rotations, hemisphere mask and optional ROI names. All models in one call must share ROI ordering, gradients and draw counts; process differing atlases separately.

## Geographic overlap

SI2 accepts original tab-separated parcel exports through `python workflows/assemble_traveling.py --first first_lh.txt first_rh.txt --second second_lh.txt second_rh.txt --output outputs/source`. Scan names contain the center identifier, such as `_C1.`, and parcel columns end in `_thickness`. Run `render_workflow` on `outputs/source/traveling_unshared/SupplementaryFigure2`.

SI4 accepts annual residence histories. MATLAB `china/generate_participant_inputs.m` exports aligned participant IDs, valid flags, residence status and birthplace/destination labels. Run `python workflows/generate_history_inputs.py --history annual.xlsx --participants participant_inputs.csv --output outputs/history`, then `python workflows/assemble_residence.py --history outputs/history/history_features.npz --participants participant_inputs.csv --output outputs/source`. Render the resulting `SupplementaryFigure4` directory with `render_workflow`. Annual workbook columns A, B, C and F contain the numeric record indicator, participant ID, year and province. Duration is the last recorded year minus the year of the first province change. Analysis uses valid rows with residence status 2.

For ED6 `reference_histograms`, pass a common reference vector or one reference per model in model-by-ROI order. Both the observed and null-map spatial correlations use Pearson; p is the inclusive absolute-tail proportion without a +1 correction. For ED10 the generic assembler accepts `full_beta` and `fold_beta` directly from `leave_one_unit_out`; add `unit_labels` in omitted-unit order.

```sh
python workflows/grid_overlap.py --input inputs/coordinates.npz --output outputs/overlap100 --size 100
```

Input `imaging_lonlat` and `culture_lonlat` are [N,2] longitude/latitude in degrees. The original Miller grid uses a window from 80掳E,60掳N to 140掳E,5掳N; occupied cells are clipped to the bundled China boundary and rendered with the Albers map projection. Outputs include cell counts, clipped polygons, SVG and PNG. Province names, geographic coordinates and weights for other maps follow `figures/README.md`. `render_maps.py` accepts optional legend entries (`label`, `color`, `size`) and inset panels (`position` in axes coordinates and a nested `panel` object).

