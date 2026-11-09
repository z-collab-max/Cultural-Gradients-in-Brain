# Study-input workflows

Run Python commands from the repository root after installing `requirements.txt`. For MATLAB adapters use `addpath('workflows/study_inputs','figures','china')`. Inputs are supplied by the caller; place resulting individual-level intermediates outside a distribution directory. Relative command arguments and MATLAB JSON paths resolve from the current working directory. Paths inside the rice-model JSON resolve from the directory containing that JSON file. Inputs to the shared China and international chains are specified in `china/SHARED_CHAIN.md` and `international/FIGURE5.md`.

China individual cortical models use the shared 28 covariates and province/city random intercepts, fitted by MATLAB REML. Python model adapters also require MATLAB. See [canonical fields and runtime requirements](CHINA_INPUTS.md). Spatial RSA/regressions and psychometric calculations retain their separate definitions. Observed effects, rotations and province shuffles must use the same named atlas order. Statistical panels save MATLAB FIG files; maps additionally save PNG files.

| Figure | Raw input and processing | Source-data output | Plot entry |
|---|---|---|---|
| Fig.1 | `generate_geography.py`: questionnaire, annual residence workbook and cultural coordinate table; `generate_mri_counts.m`: brain structure | province/city counts and grid CSVs | `figures/render_maps.py`, `workflows/grid_overlap.py` |
| Fig.2 | shared China chain; `render_fdr.m` also accepts the study effect summary | main maps, ROI means/SEM, SNR/replacement, FDR masks | shared-chain assemblers and `render_fdr.m` |
| Fig.3 | shared relocation chain; `source_data_figure3_a_agedist.m`, `source_data_figure3_e_chinamap.m`: brain structure, ComBat selection, residence history | age histograms/CDF and province cultural means | `plot_age_distribution.m`, `plot_province_culture.m`; shared bootstrap plotters |
| Fig.4 | shared behavioral/mediation model, `best_task_order.m`, Fig.4 spatial assembler | task effects, selected-order RDM/RSA/slope, spatial tests | Fig.4 specialized renderer |
| Fig.5 | `international/prepare_study_cortex.py`, `prepare_study_units.py`, `fig5_chain.py` | harmonized inputs, endpoint maps, comparisons and geographic values | `figures/render_fig5.m`, generated maps |
| ED1 | `generate_psychometrics.py`, `generate_beijing_distances.m`, `generate_brain_icc.m` | item statistics, regional distances, cortical ICC | `render_brain_icc.m`; table/distance layouts in `figures/layouts` |
| ED2 | `generate_china_trajectory_inputs.m`, `generate_china_trajectory.py`; `international/trajectory_models.py`; `assemble_center_geography.py`; threshold and simulation entries | each-PC observed/null curves, full-cohort origin/site counts, threshold comparisons | `figures/render_study_outputs.m`, generated SVG/PNG maps |
| ED3 | `china/generate_gradient_inputs.m`, `export_rebuilt_gradients.py` | nine ordered gradients and correlations | `render_gradients.m` |
| ED4 | external cortical vertex preprocessing; `run_vertex_spin.m` | vertex maps, valid masks and rotation statistics | `figures/render_vertex_map.m` |
| ED5 | `china/prepare_combat_inputs.m`, `run_combat_study_models.py`; shared mixed-model comparison entry | separately labelled protected-PC12 and protected-PC123 maps and direct-shuffle nulls | `assemble_combat_panels.py`, `figures/render_study_outputs.m` |
| ED6 | `generate_culture_effect_rsa.m` | cultural RDM, observed RSA, null and spatial comparisons | `render_culture_rsa.m` |
| ED7 | `run_area_volume.py`; `run_thickness_models.py --family pc7` | 307-ROI area/volume, 308-ROI PC3 and 19 eigenmodes, spins and province nulls | `assemble_area_volume.py`, `assemble_thickness_panels.py`, renderers |
| ED8 | `rice_culture.py`, `run_rice_models.py` | cultural correlations and four cortical rice models | `render_rice_culture.m`, `render_rice_models.m` |
| ED9 | `generate_china_resampling.m`, `export_china_resampling.py`, `international/resample_models.py` | all draws, valid-fit counts and quantiles by sampling proportion | `render_china_resampling.m`, `figures/render_panels.m` |
| ED10 | `china/run_china_deletion.m`, `international/run_loo.py` | fixed cultural residual direction, fold maps, correlations | `china/render_china_deletion.m`, `figures/render_panels.m` |
| ED11 | shared Fig.3 ROI selection, six-group aggregation and paired bootstrap | ROI/network contrasts and intervals | shared Fig.3/ED11 renderer |
| ED12 | Fig.4 behavioral and two mediation spatial branches | RSA/slope maps and their comparisons | Fig.4/ED12 specialized renderer |
| ED13 | `run_thickness_models.py --family dimensions`, `assemble_dimension_correspondence.py`; `international/dimension_correspondence.py` | six cultural-dimension China/nations comparisons and cultural/gradient correspondence | `figures/render_study_outputs.m` |
| ED14 | `run_thickness_models.py --family hof6`, `generate_hofstede_webb.m`, `assemble_thickness_panels.py --family hof6` | six-gradient fits, PC-specific spins, direct province nulls, HC2 Webb intervals | `figures/render_study_outputs.m` |
| ED15 | joint cultural PCA; `assemble_dimension_correspondence.py`; `international/dimension_correspondence.py` | joint axes/angles and PCA of six cortical dimension-effect maps | dedicated joint-PCA layouts; `figures/render_study_outputs.m` |
| SI1 | `analyses/run_word_pipeline.py` | cleaned indices, weighted scores and validation correlations | SI1 original layout |
| SI2 | traveling-participant assembler/calibration | aligned site measurements and cross-participant errors | SI2 original layout |
| SI3 | `simulation/generate_si3.m`, `export_si3_source.py` | simulation realizations and summaries | `simulation/plot_si3.m` |
| SI4 | annual residence extraction and selected-participant table | durations and geographic routes | SI4 original layout |

## Province maps and age distributions

```matlab
source_data_figure3_a_agedist('inputs/brain_data.mat','inputs/combat','inputs/residence.xlsx','outputs/figure3');
plot_age_distribution('outputs/figure3/Fig3a_AgeDist_SourceData.mat','outputs/figure3_plots');
source_data_figure3_e_chinamap('inputs/brain_data.mat','inputs/combat','inputs/residence.xlsx','outputs/figure3');
plot_province_culture('outputs/figure3/Fig3e_ChinaMap_SourceData.mat',[], 'outputs/figure3_maps');
```

Province maps use Chinese province keys. Relocators are grouped by birthplace, and their values are destination cultural scores. The geometry argument may be empty to use bundled geographic boundaries. The ComBat directory supplies `data_harmonized_thickness_demo_cul_valididx_withoutEnv.mat` with `valididx_base`.

## Area and volume

```matlab
export_birth_provinces('inputs/brain_data.mat','inputs/birth_provinces.csv');
```

```text
python workflows/study_inputs/run_area_volume.py --brain inputs/brain_data.mat --subjects inputs/birth_provinces.csv --output outputs/area_volume --draws 10000
python workflows/study_inputs/assemble_area_volume.py --source outputs/area_volume --gradients inputs/dataset_brain.mat --spins inputs/perm_idd_dk308.mat --output outputs/area_volume_panels
```

```matlab
render_area_volume('outputs/area_volume_panels','outputs/area_volume_plots');
```

The brain structure contains named area/volume ROI columns, demographic/genetic/environmental fields listed in the generator, cultural PCs, recruitment scores and center labels. The keyed CSV contains `subject_id,province`. The point-estimate design and null design are saved separately; `null_model_effect` is the unshuffled estimate under the null design. Rotations are a 308-by-B one-based integer matrix `perm_id`; 307 retained rows still index all 308 source effects. PC2 changes orientation at assembly. Its displayed spin regression omits eigenmode 2, while its province regression includes all six gradients.

## International sensitivity analyses

All entries accept the same raw input folder as Fig.5. Optional `--harmonized` reuses saved arrays from that chain.

```text
python international/threshold_maps.py --inputs inputs/international --output outputs/thresholds
python international/resample_models.py --inputs inputs/international --output outputs/resampling --draws 100
python international/run_loo.py --inputs inputs/international --output outputs/loo
python international/dimension_maps.py --inputs inputs/international --dimensions inputs/cultural_dimensions.csv --components 0 0 --output outputs/dimensions
```

The dimension table contains `unit` plus named numeric cultural dimensions. Component numbers are zero-based and must be chosen explicitly. Threshold comparison defaults to endpoint; `--reference 99.9` requests a separate threshold reference. Individual resampling uses all final-model covariates at proportions 10–90%. Whole-unit resampling uses demographic within/between columns at proportions 50–90%; the full reference is refitted with that same reduced design. Each failed rank-deficient fit is saved explicitly with valid-draw counts. LOO fixes full-sample cultural residual directions and rebuilds artifacts within each fold.

```matlab
render_panels('outputs/loo/PC1/ExtendedDataFigure10/panels.mat','outputs/plots');
run_china_deletion('inputs/china','outputs/china_chain','outputs/china_deletion',true);
```

## Other raw inputs

- `generate_beijing_distances(culture_file,extra_file,out,draws,seed)` expects the original cultural questionnaire structures and supplementary province records. Regional distances use the reference-province variance; bootstrap variances use paired Beijing/current-region samples.
- `generate_brain_icc(brain_file,out)` computes cortical adjustment and province ICC. `generate_psychometrics.py --help` lists questionnaire and response-scale options.
- `run_vertex_spin(pkg,enigma,out)` expects prepared vertex-effect and biological-gradient MAT files, fsaverage5 masks, ENIGMA sphere surfaces and `SurfStatReadSurf1`. Vertex smoothing and individual cortical reconstruction precede this entry point.
- `generate_culture_effect_rsa(config,out)` requires JSON fields `brain_file,raw_file,effect_file,dataset_file,spins_file`; optional `seed,draws,rsa_statistic,surface_statistic` select the run. Observed and shuffled RDM correlations both use `rsa_statistic` (default Spearman); observed and null-map reference correlations both use `surface_statistic` (default Pearson). Culture-effect scatter spins use Pearson. Province shuffles are unrestricted.
- `generate_china_resampling(config,out)` requires the study raw-result and brain inputs named in its config reader. It defaults to 100 replicates at each of nine proportions. Raw individual arrays are required to refit subsets.
- Rice adapters accept a province table and original indicator SAV, or the aligned province fields specified by `--help`. `run_rice_models.py` requires ID-keyed subject linkage, raw brain fields, gradients and rotations. Its JSON contains `brain`, `subjects`, `indicators`, `roi_order`, `gradients`, `spins`; all six paths are relative to the JSON directory. `subjects` is a CSV with `subject_id,province`; `indicators` is a CSV with `Province_CN,PercentPaddy,EnvironmentalRiceSuitability`. The ROI MAT contains `Talhelm_Gradient_Results/thickness/brain_var` and `Talhelm_Gradient_Results/area/brain_var` in paired atlas order; `perm_id` contains one-based ROI-by-rotation indices.
- `assemble_alternative_gradients.py` requires the **alternative six-dimension culture** model's `beta,null_beta,cluster_influence,roi_names`, not an unrelated seven-dimension model. Its spatial NPZ must supply exactly matching `roi_names`, named gradients, and zero-based rotation-by-ROI indices. Run the shared China chain with the appropriate keyed cultural PC table first.

To inspect arguments use `python <entry.py> --help`. See `workflows/SPECIALIZED_INPUTS.md` for the remaining named-array schemas and original-layout renderers. Input availability and statistical equivalence must be checked separately from successful figure export.
## Additional dimension-omission models

A completed shared China run supplies full-sample ComBat data. Each cultural dimension and each leave-one-dimension-out pair then receives its own artifact endpoint and 28-covariate, province/city-intercept mixed model. ComBat is fixed across these dimension comparisons.

```matlab
run_dimension_branches('inputs/china','outputs/china_chain','inputs/dimensions.csv','outputs/china_dimensions',10000,42);
```

The dimension CSV has `province` as its first column and finite numeric cultural dimensions in the remaining columns. It covers all birth and center provinces. Leave-one-dimension-out PCs are matched one-to-one to the two full cultural PCs by maximal absolute correlation, with orientation saved explicitly. Output `dimension_models.mat` includes full/fold maps, dimension-effect maps, ROI names, component matching and signs.

```text
python workflows/assemble_dimension_branches.py --models outputs/china_dimensions/dimension_models.mat --spatial inputs/spatial.npz --components 0 1 --output outputs/china_dimension_source
```

The spatial NPZ contains identical `roi_names` and zero-based **rotation-by-ROI** `spin_indices`. Effect-PCA components are selected explicitly by the two zero-based `--components` values. Render both output `panels.mat` files using `figures/render_panels.m`.


## Thickness, cultural dimensions and Webb intervals

The study mixed-model adapters use canonical `brain_data` fields and the `Culdata` province table. The 28 covariates and province/city group fields are specified in [CHINA_INPUTS.md](CHINA_INPUTS.md). Named atlas ROI fields are selected in the bundled order. Optional `--indices` accepts `indices_zero_based` [province,draw] into `Culdata`; otherwise shuffles are generated with the saved seed. No outer correlation screening is applied.

```text
python workflows/study_inputs/run_thickness_models.py --brain inputs/brain_data.mat --culture inputs/Culdata.mat --family pc7 --output outputs/thickness_pc7 --draws 10000
python workflows/study_inputs/run_thickness_models.py --brain inputs/brain_data.mat --culture inputs/Culdata.mat --family hof6 --output outputs/thickness_hof6 --draws 10000
python workflows/study_inputs/run_thickness_models.py --brain inputs/brain_data.mat --culture inputs/Culdata.mat --family dimensions --output outputs/thickness_dimensions --draws 10000
```

```matlab
addpath('workflows/study_inputs','china','figures');
generate_hofstede_webb('inputs/brain_data.mat','inputs/dataset_brain.mat','outputs/hofstede_webb.mat');
```

```text
python workflows/study_inputs/assemble_thickness_panels.py --source outputs/thickness_pc7 --family pc7 --gradients inputs/dataset_brain.mat --output outputs/pc7_panels
python workflows/study_inputs/assemble_thickness_panels.py --source outputs/thickness_hof6 --family hof6 --gradients inputs/dataset_brain.mat --spins inputs/perm_idd_dk308.mat --webb outputs/hofstede_webb.mat --output outputs/hof6_panels
python workflows/study_inputs/assemble_dimension_correspondence.py --dimensions outputs/thickness_dimensions --hofstede outputs/thickness_hof6 --gradients inputs/dataset_brain.mat --spins inputs/perm_idd_dk308.mat --output outputs/china_correspondence
```

The PC3 surface uses the REML mixed-model t statistic and saves its degrees of freedom. Webb uses 10,000 balanced antithetic draws with MATLAB Twister seed 42 and province HC2 influence scores from the fitted mixed-model covariance. The MAT bridge retains the exact weights. PC2 spatial rotations omit eigenmode 2; province tests retain all six gradients.

## ComBat study branches

```matlab
prepare_combat_inputs('inputs/brain_data.mat','inputs/Culdata.mat','outputs/combat_inputs.mat');
```

```text
python workflows/study_inputs/run_combat_study_models.py --inputs outputs/combat_inputs.mat --protected-pcs 2 --null-protected-pcs 2 --draws 10000 --output outputs/combat_PC12
python workflows/study_inputs/run_combat_study_models.py --inputs outputs/combat_inputs.mat --protected-pcs 3 --null-protected-pcs 3 --draws 10000 --output outputs/combat_PC123
python workflows/study_inputs/assemble_combat_panels.py --source outputs/combat_PC12 --protected-pcs 2 --gradients inputs/dataset_brain.mat --spins inputs/perm_idd_dk308.mat --output outputs/combat_PC12_panels
python workflows/study_inputs/assemble_combat_panels.py --source outputs/combat_PC123 --protected-pcs 3 --gradients inputs/dataset_brain.mat --spins inputs/perm_idd_dk308.mat --output outputs/combat_PC123_panels
```

Observed and null fits use the same protected-PC count. The optional `--null-protected-pcs` value defaults to `--protected-pcs`; mismatched values are rejected. The no-culture, PC12-protected and PC123-protected matrices are all generated in the preparation step. The shared mixed-model version remains `china/run_combat_comparison.m`.

## Empirical trajectories

```matlab
generate_china_trajectory_inputs('inputs/brain_data.mat','inputs/province.xls','inputs/center_shuffles.mat','outputs/china_trajectory_inputs');
```

The province workbook uses `Sheet2` with `pc7_1` and `pc7_2`. The shuffle MAT supplies `center_ids` and `culture_labels` [center,at least 10000 draws], with one-based province-table row indices. The generator saves the individual selection and 29 residual matrices per PC. Removal counts are 0 and 2–29. The inner |r| < .01 filter is applied before each ROI-specific artifact PCA.

```text
python workflows/study_inputs/generate_china_trajectory.py --brain-file inputs/brain_data.mat --raw-file outputs/china_trajectory_inputs/raw_res_birth_all_expcs.mat --culture-file inputs/Culdata.mat --index-file inputs/province_indices.npz --target-file inputs/dataset_brain.mat --draws 2000 --output-dir outputs/china_trajectory
python international/trajectory_models.py --inputs inputs/international --output outputs/international_trajectory --draws 1000 --verify
python international/trajectory_extra.py --inputs inputs/international --output outputs/extra_diagnostic
```

International trajectories refit each count from zero through endpoint for every shuffled cultural score. Checkpoints contain input hashes and seeds. `--verify` compares cross-product fits against the direct protected fit. `trajectory_extra.py` projects out the full unprotected score space, saves effect norms and ranks, and returns an undefined spatial correlation when the map norm is below 1e-10. Optional `--draws` repeats that diagnostic on shuffles.

## International dimension correspondence

```text
python international/dimension_correspondence.py --inputs inputs/international --dimensions inputs/dimensions.csv --china-maps outputs/thickness_dimensions/models.npz --gradients inputs/dataset_brain.mat --output outputs/international_dimensions --draws 10000
```

`dimensions.csv` has `unit,PDI,IDV,MAS,UAI,LTO,IVR`. Country keys match Fig.5 inputs. China ROI names are matched by atlas label after removing the thickness suffix. Each dimension has its own protected endpoint and direct-shuffle null. Six comparisons, the gradient correspondence matrix, effect-PCA biplot and effect/culture-PC comparisons are assembled. Each culture PC is compared with effect PC1. `--dimension PDI` runs one dimension independently; after all six runs, `--assemble-only` builds the combined outputs from their completed checkpoints.

## Origin and scanning-center maps

```text
python workflows/study_inputs/prepare_geography_participants.py --nations inputs/nations_cortex.mat --appended inputs/usa_append.mat --output outputs/nations_geography_participants.csv
python workflows/study_inputs/assemble_center_geography.py --brain inputs/brain_data.mat --nations-participants outputs/nations_geography_participants.csv --nations-coordinates inputs/capital_coordinates.csv --usa-participants inputs/usa_cohort.csv --output outputs/geography
```

The original cortical identifiers determine country and dataset after dataset reassignment and duplicate-source removal. Geographic summaries use the full input cohort. The country-coordinate CSV has `unit,latitude,longitude` (or `nation,cap_lat,cap_lon`). USA rows contain `state,dataset,latitude,longitude`. China origin coordinates are assigned to the bundled province polygons and summarized by scanner; centers use the modal origin province location. Only aggregate counts and coordinates enter the map recipe.

`generate_geography.py` normalizes the four municipality city labels and combines Panzhou with Liupanshui for city aggregation. Coordinates must be joined by city key, never by workbook row position.

## Render completed runs

```matlab
addpath('figures');
render_study_outputs('outputs/hof6_panels','outputs/hof6_figures');
```

The renderer walks `panels.mat` and named `*surface.mat` inputs and saves FIG plus PNG. Run it once for each analysis output directory. All statistical axes use Arial 8 and line width 0.85.
