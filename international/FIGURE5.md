# Figure 5 from individual cortical measurements

Run from the repository root after installing `requirements.txt`.

```sh
python international/make_fig5_inputs.py --output outputs/fig5_inputs
python international/fig5_chain.py --inputs outputs/fig5_inputs --output outputs/fig5 --permutations 4 --webb-draws 4
```

The example creates synthetic participants in 26 countries and 18 US states and intentionally reorders input rows to exercise the joins. The small draw counts and synthetic ROI permutations demonstrate the interface. For inference, supply valid atlas rotations and set the draws explicitly. Defaults are 10,000 unrestricted unit shuffles and 1,000 Webb refits. Use a new output directory for each run.

Optional `--models PC1` runs one model. Permutation blocks retain input hashes, seeds and draw counts. Use `--resume` to reuse completed blocks with identical settings; observed fits and Webb resampling are recomputed. The cross-product endpoint implementation checks its first three shuffled fits against the direct protected fit.

Optional `--permutation-indices inputs/outer_permutations.npz` supplies fixed outer permutations. Each selected model has an integer array named `PC1`, `PC2` or `USA_TL` of shape `[permutations, units]`, with a zero-based permutation in each row. The matching string arrays `PC1_unit_names`, `PC2_unit_names` and `USA_TL_unit_names` identify the unit order. Supplied indices take precedence over the random seed. Cache validation includes the null model and supplied indices.

```matlab
addpath('figures');
render_fig5('outputs/fig5/source','outputs/fig5/figures',false);
```

Use `true` as the final argument for the 308-region atlas to also draw the cortical t maps. Statistical panels save editable `.fig` files. The Python run writes each geographic map as SVG and PNG.

## Input files

The input directory contains `Nations` and `USA` subdirectories. Participant IDs may repeat across sites, so supply a unique row key for each scan (for example, include original row, site and participant ID). Use the included sample, with complete age, modeled sex and TIV values, including flags indicating which age and sex values were observed. Cortical thickness is in its original measurement units.

| File | Contents |
|---|---|
| `<cohort>/participants.csv` | `row_key`, `unit`, `batch`, `age`, `female`, `tiv`, `age_observed`, `sex_observed`. The two flags use 0/1. |
| `<cohort>/cortex.csv` | `row_key` and one original thickness column per ROI. These column names define the ROI order. Rows are matched by key. |
| `<cohort>/units.csv` | One row per named `unit`, with culture and macro fields below. USA also supplies `abbreviation`. |
| `<cohort>/pool.npz` | `pool` [draw × unit], zero-based permutations for artifact estimation, and `unit_names` in first-occurrence order of the modeled units. |
| `targets.csv` | Unique `roi` with named external spatial reference columns `eigenmode2`, `SA`, `China_PC1`, `China_PC2`, `China_TL`. Supply the intended sign of each target. Rows are aligned by ROI name. |
| `spins.npz` | `indices` [ROI × spin], zero-based spatial rotation indices; `roi_names` must match the model's ordered ROIs. |

Nations culture columns are `PC1`, `PC2` in country score units. Macro fields follow `NATIONS_BLOCKS` in `core/model_spec.py`: `HDI2020`, `HDI1990`, `PPPrw`, `tmp`, `frs`, `wet`, `freshwater`, `arable_land`, `food_production`, `abslatitude`, `PM2_5`, `forest_cover`. USA uses raw `USA_TL` and the eight `USA_MACROS`: `humidity`, `disaster_rate_per_people`, `land_forested`, `population_density`, `precipitation`, `air_quality_index`, `PPP_2017`, `highest_elevation_ft`.

## Original study-table adapter

For the study's native MAT/CSV schemas, the alternative first step is:

```sh
python international/prepare_study_cortex.py --config inputs/cortex_paths.json --output outputs/harmonized
python international/prepare_study_units.py --config inputs/unit_paths.json --harmonized outputs/harmonized --output inputs/figure5
python international/fig5_chain.py --inputs inputs/figure5 --harmonized outputs/harmonized --output outputs/fig5_study --permutations 10000 --webb-draws 1000
```

The JSON names nine paths, relative to the JSON file: `nations_cortical_table`, `nations_appended_usa_table`, `usa_cortical_table`, `usa_batch3_cortical_table`, `nations_demographics`, `usa_demographics`, `demographic_changes`, `demographic_decisions`, `pnas_identity_resolution`. The first two are the Nations `nationsdata_all` MATLAB structure and US `roi_names/state/thickness/subid` MATLAB table. The USA cortical tables have `subid`, `state`, `dataset` and named ROI columns. Demographic manifests provide original row/participant or subject key, cultural unit, scanner batch, demographics and TIV; the three JSON inputs provide the original correction records, exclusions and participant-identity decisions. This adapter applies the study-specific site transformations, appends the three US state datasets to Nations, merges the additional USA cohort, and verifies the original Nations row/identity pairs and USA subject/site matches. With `--harmonized`, the cortical CSV and participant CSV in the generic schema are superseded by the adapter outputs. The unit tables, pool, reference maps and spatial rotations remain required.

The unit adapter accepts `nations_macro` (NPZ with `g`, `macro`, `macro_labels`, `C`, `pool`, `targets`), `nations_pc1` (NPZ with `g`, `X`, `labels`, `C`, `pool`, `targets`, including the `PPPrw` column), `usa_macro` (`C`, `macro_group`, `labels`), and `usa_model` (`g`, `C`, `pool`, `targets`). `nation_culture` is a MAT file containing `matched_nations` and the local country PCA `score_world`. `state_culture` is the `usa_current` workbook with `State`, `Sname`, `Cultural_Tightness`. Supply the declared reference `roi_order` as one parcel per line and `spin_indices` as a zero-based NPY matrix. All paths in `unit_paths.json` are relative to that file. The adapter checks individual cultural-unit coding, constancy of country-level macros, the standardized cultural scores, common eigenmode values and parcel ordering before writing model inputs.

## Calculation and outputs

ComBat is applied within each cultural unit. The design uses the cohort-specific individual-within-unit and between-unit demographic terms. Nations standardizes the twelve macro measures across countries and retains all PCA components within each named macro block. USA standardizes the eight state macro measures at the participant level after mapping states. The observed nuisance designs have 29 and 21 columns, respectively. Artifact construction retains its internal `|r| < .01` screening against the current candidate culture. Endpoint uses `number of units − 1` PCs. Each outer draw shuffles unit culture scores without an outer correlation filter. The null nuisance design uses the same covariates as the observed model; the shuffled score is the predictor of interest. Observed effects use the original design. Webb resamples refit endpoint artifacts and protected REML. The ROI effect map and both spatial comparisons are evaluated for each refit.

Comparison 1 uses eigenmode 2 for Nations PC1 and USA TL, and the S–A axis for Nations PC2. Comparison 2 uses the corresponding Chinese cultural effect. Spin rotates the reference, keeping the observed international effect map fixed. Spin and permutation probabilities use the absolute tail without a +1 correction. Comparison files retain `p_permutation` and `permutation_null`, with `null_method` identifying the unit-shuffle test. The CI uses the 2.5 and 97.5 percentiles of the Webb correlation draws. Displayed cortical maps use t values; scatterplots use standardized coefficients. Geographic source data preserve raw and model scores separately. Nations PC1 is sign-reversed only for geographic display. The original palettes use −4 to 4 for countries and 20 to 70 for state tightness.

`source/` contains the brain, comparison and geography MAT/NPZ files for each model plus map previews. `local/` contains the harmonized individual data, model designs, cleaned values and resampling details. Keep individual-level outputs in the analysis environment.
