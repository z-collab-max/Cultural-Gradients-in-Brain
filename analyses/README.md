# Additional methods

Run the synthetic numerical examples from the repository root:

```text
python analyses/run_methods_example.py
```

For the raw daily word-index workflow, supply the workbook and province table:

```text
python analyses/run_word_pipeline.py --raw-input inputs/daily_words.xlsx --province-input inputs/province_scores.csv --output outputs/word_analysis
```

The workbook sheet `Sheet1` needs `省/直辖市`, `关键词`, `日期`, and `出现关键词的微博指数`. The province CSV needs `names` and `z_tightness_control`. The function applies the specified daily-index replacements, aggregates the 20 keywords, aligns 31 provinces and computes the variance-weighted PCA scores. The retained-component polarity settings are `[1,1,1,1]` and `[1,-1]`; a mismatch in retained counts raises an error. Input files are not distributed.

They save MAT/NPZ source arrays under `outputs/additional_methods`. For Python use in a custom script, add `analyses` to `sys.path` before importing the modules. For MATLAB, use `addpath('china'); addpath('analyses')`. See [the figure index](FIGURE_COVERAGE.md) for inputs and coverage boundaries.

The exported psychometric, weighted-word PCA, traveling calibration, map-regression and LOO protection kernels retain the corresponding original numerical calculations. Other modules provide explicit input/output interfaces for related calculations.

## Array contracts

| Function | Input contract | Output |
|---|---|---|
| `anova_icc` | Response vector and group IDs | F, p, ICC1, ICC2, sample and group counts |
| `rwg` | DataFrame with `from_province`, item columns, response-scale maximum | Mean rwg and per-group records |
| `cronbach_alpha` | Complete observations × items matrix | Alpha |
| `cultural_pca` | Units × cultural dimensions | Scores, loadings, explained variance and scaling |
| `joint_pca_axis_angles` | Joint-to-local score-correlation matrices for China and world | Acute angles within the displayed first-two-PC plane |
| `weighted_score` | Normalized word frequencies and one polarity per retained PC | Variance-weighted score and PCA metadata |
| `commonality` | Stayer/relocator difference map and three aligned component maps | Unique/shared/incremental R² components, coefficients, fitted map |
| `shared_roi_bootstrap` | Positive and negative individual ROI summaries, group IDs | Minimum-group-size sampling with replacement; identical indices for both summaries |
| `spin_test` | Two ROI maps, rotations × ROIs integer index array | Correlation, null distribution, absolute-tail two-sided p |
| `partial_map_coefficients` | ROI effect vector and ROI × gradients matrix | Simultaneous standardized gradient coefficients |
| `harmonize` | Individuals × ROIs cortex, batch IDs, numeric protected covariates without intercept | ComBat-adjusted matrix |
| `remove_artifacts` | Individual cortex; individuals × shuffles scores; shuffles × ROIs effects; candidate-relative pool correlations; component count | Cleaned cortex and explained variance |
| `wild_cluster_draws` | Fitted values, residuals, cluster IDs, statistic callback | Webb or Rademacher bootstrap statistics; caller defines interval/studentization |
| `snr_curve` | Matched signal/residual matrices, group IDs, fitted-map callback, positive signal/noise mean-square power ratios (variance ratios for centered inputs) | Fitted maps across levels and repetitions; residual rows are resampled within group; noise amplitude scales as the inverse square root of SNR |
| `fit_predict` | Training reference/site ROI vectors, held-out reference vector, model name | Held-out calibrated measurements |
| `leave_one_unit_out` | Dictionary with `X`, `Y`, `g`, `C`, `pool`, `targets` | Full effect map, fold maps, map correlations, artifact PC counts |

`g` and Python spin/permutation arrays are zero-based. MATLAB mesh faces are one-based. Spin correspondences must come from valid atlas/spherical geometry; arbitrary row shuffles are not spatial rotations. `grid_overlap` accepts already-projected coordinates and a cell size in the same unit; coordinate projection and boundary clipping are external inputs.

For raw China inputs, use `prepare_china_table` before fitting. `null_effect_coefficients` uses the common 28-covariate/two-intercept specification. `run_sensitivity` adds named covariate sets to that specification. `refit_sample` accepts distinct participant indices for deletion/subsample analyses. All outer culture shuffles continue to use `permutation_china`.

For traveling-participant calibration, valid models are `Additive`, `Multiplicative`, `ROI-specific additive`, `Affine`, `Polynomial2`, `Polynomial3`, and a radial-basis model. Keep the participant used to estimate calibration separate from the participant used to evaluate it.

For the LOO calculation, covariate scaling and ComBat are held fixed; cultural residual directions are estimated on the full sample. Each fold rebuilds its artifact basis, updates weights, and refits REML. The number of retained artifact PCs is the fold's number of units minus one.

The numerical kernels do not identify an empirical cohort or silently retrieve missing inputs. Supply arrays with complete, aligned rows and columns. Large bootstrap/permutation jobs should be checked on a small subset before a full run.
