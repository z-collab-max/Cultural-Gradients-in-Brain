# Shared China cortical preprocessing and model fits

From the repository root in MATLAB:

```matlab
addpath('china');
base = fullfile(pwd,'outputs','china_shared');
make_china_chain_inputs(fullfile(base,'inputs'));
run_china_chain(fullfile(base,'inputs'),fullfile(base,'results'),3,1200,42);
```

The example starts with synthetic participant-level regional measurements. The three outer draws demonstrate data flow. For inference, set the outer and inner draw counts explicitly (both default to 10,000). Each outer draw repeats harmonization, artifact estimation, the main model and the relocation models. Choose a new results directory for each run.

| File | Fields |
|---|---|
| `participants.csv` | `participant_id`, `valid` (0/1 inclusion mask), `province_birth`, `city_birth`, `province_destination`, `relocator` (0/1), `scanner` (ComBat batch), `center_province` (center-to-cultural-unit mapping for artifact estimation) |
| `covariates.csv` | `participant_id` and the 28 numeric covariates in `china_covariates.m` |
| `cortex.mat` | `participant_id`, numeric `cortex` [participants × ROIs] in raw measurement units, ordered `roi_names` |
| `culture.csv` | One row per province with `province`, `pc1`, `pc2`; the two cultural scores must already have their intended orientation |

IDs are matched before any fit. Scanner IDs, center-to-province labels and birthplace IDs have distinct roles. Supply an eligible participant mask, finite analysis values and at least two observations per scanner. The protected and nuisance designs must have full rank. Destination and center provinces must be represented among eligible birthplace provinces. The example shuffles raw cortical and covariate row order independently of the participant table.

An original `brain_data` structure in `cortex.mat` is also accepted. It must contain `MRIID` and one numeric measurement vector per parcel. Supply `roi_names.txt` listing the exact field names in atlas order, one per line. Both cortical input formats use the same participant-ID match.

1. Align the eligible rows. Apply parametric ComBat to raw cortical measurements by scanner, protecting the 28 covariates, both birth scores and both destination scores.
2. Generate center-level shuffles of cultural-unit labels. For each PC separately, retain candidate predictors with absolute individual-level correlation below 0.01 against that draw's target. Estimate their ROI coefficients by OLS adjusted for the covariates. Center the candidate contributions, weight by the number of individuals per center, and remove the endpoint components separately for each PC and ROI. The endpoint is the number of center-level units minus one. The code requires enough accepted candidates and full endpoint rank.
3. Fit the main cultural association and relocation models using `china_formula.m`, with the 28 covariates and random intercepts for birthplace province and city within province. The relocation fits include birthplace culture, relocation status, their interaction, destination change and separate group fits.
4. For every outer draw, shuffle paired PC rows across birthplace provinces, map the same shuffled rows to birthplace and destination, and repeat the preceding steps. The shuffle is unrestricted. All candidates share the inner center orders. The ROI permutation p value uses `mean(abs(null) >= abs(observed))`, without a +1 correction.

`main_maps.mat` stores PC-by-ROI point estimates, model SE/p, permutation p, `null_beta` [outer draw × PC × ROI], and observed HC2 province influence for Webb intervals. `commonality.mat` stores the five map families for the Figure 3 decomposition, whose target is mover minus stayer. `relocation_results.mat` contains the fitted and shuffled relocation coefficients and group maps. `preprocessing.mat` contains the aligned IDs, harmonized measurements, PC-specific clean measurements, applied orders, inner selection counts and observed inner selection. Keep individual-level outputs in the analysis environment.

To build the gradient panels and commonality panels, provide a spatial NPZ containing ordered `roi_names`, `gradients` [ROI × gradient] and zero-based `spin_indices` [rotation × ROI]. Optional fields are `gradient_labels` and `left_hemisphere`. Supply spatial rotations appropriate to the atlas.

```sh
python workflows/assemble_china_chain.py --results outputs/china_shared/results --spatial inputs/spatial.npz --output outputs/china_shared/source --repetitions 10000
```

```matlab
addpath('figures');
render_panels(fullfile(base,'source','Figure2','panels.mat'),fullfile(base,'figures'));
render_workflow(fullfile(base,'source','commonality','Figure3'),fullfile(base,'figures'));
```

## Figure 2 provincial scatter and recovery experiments

New shared runs write `fig2_individual_inputs.mat` alongside `main_maps.mat` in the run directory. This local file contains raw thickness (mm), ordered participant/ROI IDs, scanner, birthplace province and city (`province`, `city`), birth PC scores and the aligned 28-covariate design. Regenerate older outputs without `city`; no substitute grouping is inferred. Run the same shared chain with the scatter cohort mask if the scatter and recovery cohorts differ. Pass its results using `--scatter-results`.

```sh
python workflows/assemble_fig2_details.py --results outputs/china_shared/results --output outputs/fig2_details --targets ROI1 ROI2 ROI3 ROI4 --draws 30 --test-orders 60
```

The four ROI arguments in this command are the synthetic example labels. For the study's atlas omit `--targets` to select left superior frontal part 8, left lateral occipital part 7, left superior frontal part 9 and right inferior parietal part 4 in PC1 positive/negative and PC2 positive/negative order. Omit `--draws` and `--test-orders` for 1,000 draws in each of two independent batches and 2,000 random-order tests per PC. The default scatter annotates the shared model's beta and p. Use `--annotations` with a JSON array of `panel`, `beta`, `p` entries in `pc1_pos`, `pc1_neg`, `pc2_pos`, `pc2_neg` order to display supplied annotations.

The scatter removes scanner means and standardizes residual thickness using sample SD. Each selected ROI is fitted with MATLAB REML using its cultural PC, all 28 covariates, and province/city random intercepts from `china_formula.m`. It subtracts only the fitted fixed covariate contributions, retaining the intercept, cultural effect and group variation, then restores the raw parcel mean and sample SD in mm. It averages within birthplace province and plots vertical SEM using population SD divided by the square root of the province sample count. The PC2 scatter uses the inverted display axis. Annotations come from the shared main-model results (or the explicitly supplied annotation file). Python plotting assembly requires MATLAB on PATH or `MATLAB_EXECUTABLE`. Counts and point coordinates are saved in `roi/Figure2/province_points.csv`.

The recovery experiment uses scanner-centered raw thickness in micrometers. Baseline fits and every bootstrap refit use MATLAB REML with the shared 28 covariates, cultural predictor, and province/city random intercepts. It centers conditional residuals within birthplace city, scales each parcel's residual pool to its fitted residual variance, and resamples participant rows within city with common indices across parcels. Independent province and city random intercepts are generated with their fitted marginal variances and cross-parcel correlations estimated from their BLUPs (50% identity shrinkage). Each city needs at least two observations, and the cohort needs multiple cities within provinces. Coefficient errors from one calibration batch are scaled against the shared main-model SE, and those factors are applied to an independent validation batch. SNR uses signal variance divided by the sum of residual, province and city noise variances. The two PCs share an SNR grid containing both actual baselines. Fixed-design recovery curves add target coefficient maps to these refitted error maps. Replacement ranks the standardized target by absolute t and replaces selected values with its signed median; it uses 100 random comparison orders, a lower-tail test without a +1 correction at nine deciles per PC, and Holm adjustment across all 18 tests. `draws/` saves both error batches, all three fitted variance components, seeds, main-model targets and calibration scales. MATLAB is required for these fits; full default runs refit both variance components for every draw and parcel.

```matlab
for branch={'roi','snr','replacement'}
    render_workflow(fullfile('outputs','fig2_details',branch{1},'Figure2'),fullfile('outputs','fig2_figures'));
end
```

## Figure 3 and Extended Data Figure 11 ROI analyses

Supply `roi_spatial.mat` with `roi_names` in the model's ROI order, `gradients` [ROI × 2] containing negative eigenmode 2 and the S–A axis, and `pc_sign` [1 × 2] containing the multipliers from model to display orientation. Use `[1,-1]` when both cultural scores retain the original PCA signs, or `[1,1]` for scores already in display orientation. The same signs apply to the birthplace model coefficients and both cultural-score tables; the gradient columns must already be oriented.

```matlab
addpath('china'); addpath('figures');
run_fig3_roi_chain(fullfile(base,'inputs'),fullfile(base,'results'), ...
    'roi_spatial.mat',fullfile(base,'roi_analysis'),10000);
render_fig3_roi_chain(fullfile(base,'roi_analysis'),fullfile(base,'roi_figures'),true);
```

`true` requests surface figures and requires real atlas ROI names. For synthetic ROI names use `false` and a synthetic pair of spatial references that selects at least one ROI per sign. Simulated rotations or gradients only exercise the interface.

The positive set requires gradient > 66th percentile and birthplace β > 90th percentile. The negative set requires gradient < 33rd percentile and birthplace β < 10th percentile. Percentiles are taken separately across ROIs for each PC. The shared run supplies the birthplace term from the relocation mixed model. Harmonized thickness is residualized against the shared covariate design and standardized per ROI across participants before taking the arithmetic mean within each signed ROI set.

Classify both birthplace and destination using the participant-weighted median of birthplace scores. PC1 ties are high and PC2 ties are low. HH and LL are stayers remaining high and low, respectively. HL and LH are movers who cross the threshold in opposite directions. MHH and MLL are movers who stay on the same side of the threshold.

Fig.3j bootstraps HH and LL with replacement, drawing the smaller of the HL and LH counts from each. Positive and negative ROI means share participant indices. The contrast is positive minus negative pooled Cohen's d, with an inclusive two-sided equal-tail p value. ED11a instead uses each participant's positive-minus-negative network mean and the smaller strict bootstrap tail with a `1/B` floor. Fig.3k uses the rounded mean HL/LH count and the original equally weighted variance denominator for each ROI's effect size. ED11b shows the four observed group means; the 95% percentile intervals resample the HH and LL means separately using each group's actual size. The forest panels use separate fixed seeds from the two network-level bootstraps.

`roi_source.mat` contains the ROI masks, thresholds, labels, per-ROI plotting arrays and Fig.3j draw vectors. `network_source.mat` contains the ED11a draws, group summaries and panel recipes. `individual_roi_means.mat` contains aligned IDs, group codes (HH, LL, HL, LH, MHH, MLL in order; zero for unassigned) and participant-level positive/negative means. Keep the last file with local analysis data.

## Figure 4 and Extended Data Figure 12 spatial tests

Use `run_china_example` to generate a behavior and mediation result with enough regional measurements for a full-rank six-gradient design. For example, `run_china_example(fullfile(pwd,'outputs','china_hierarchy'),0,10)` fits ten synthetic ROIs with point-estimate indirect effects. A bootstrap fit uses indirect effects divided by their bootstrap SE. Generate the spatial inputs for the same ordered ROI rows.

```matlab
addpath('china'); addpath('figures');
run_fig4_spatial_chain('outputs/china_hierarchy/source_data/hierarchy_inputs.mat', ...
    'inputs/hierarchy_spatial.mat','outputs/hierarchy_spatial');
render_fig4_spatial('outputs/hierarchy_spatial','outputs/hierarchy_figures',false);
```

Supply `roi_names`, `gradients` [R,6], `gradient_names`, `spin_indices` [R,S] with **one-based MATLAB ROI indices**, and logical `left_hemisphere` [R,1]. Name the gradients `SA_axis`, `developmental_expansion`, `evolutionary_expansion`, `intersubject_variability`, `cognitive_PC1` and `eigenmode2`, using the intended atlas orientation. Spins are externally computed spatial rotations. Model ROI identifiers are included in `hierarchy_inputs.mat`; for an earlier file, supply its independently recorded ROI order as `model_roi_names` in the spatial MAT. Optional `example_names` is a three-by-two cell array of ROI names for behavior, PC1 and PC2. Otherwise the two extremes of the task-slope map are used. Set the render flag to `true` for ROI names matching the bundled surface atlas.

Each branch selects a single best-fitting order. Fig.4 fits six standardized gradients jointly to the RSA map, rotates the RSA response for the spin reference and refits every retained task-order map for the task-order reference. It also saves the zero-order Pearson correlation and both null distributions. Behavior excludes the 24 within-category permutations of the selected category sequence. The mediation branches use all 5,040 orders. All coefficient and zero-order spin/task-order tests use inclusive absolute exceedances, `mean(abs(null) >= abs(observed))`, without a +1 correction.

ED12 saves the seven-task regression slope and the example task values separately from the slope used in the RSA–slope correlation. For behavior this comparison uses the task slope. For the two mediation branches it uses the slope fitted to the 21 signed lower-triangle RDM differences. The spatial test rotates the slope, keeps the RSA map fixed and counts absolute exceedances. `fig4_source.mat` and `ed12_source.mat` in each branch folder retain the test statistics, nulls, atlas labels and orders. `panels.mat` is the input for the editable statistical plots.
