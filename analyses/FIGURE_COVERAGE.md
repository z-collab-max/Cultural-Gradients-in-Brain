# Figure-to-code map

This table maps figure analyses to reusable calculation functions. It is a method-level index, not a claim that every published panel can be reconstructed from the distributed inputs. Only simulation references and synthetic inputs are distributed. Empirical row selection, variable alignment, atlas geometry, spin indices and approved aggregate values must be supplied by the caller where indicated.

Paths are relative to the repository root. China individual cortical models, including the study-input adapters, use `china/china_formula.m` with 28 covariates and province/city random intercepts. See [China input requirements](../workflows/study_inputs/CHINA_INPUTS.md). Cortex matrices can contain ROI or vertex columns. International examples use the endpoint. Outer permutations use direct cultural-score shuffles; inner artifact construction retains its correlation screen.

| Figure | Calculation entry points | Required inputs beyond the synthetic examples |
|---|---|---|
| Figure 1 | `analyses/group_summaries.py`; `workflows/grid_overlap.py`; `figures/render_maps.py` | Residence longitude/latitude, province IDs; grid projection, clipping and boundary geometry are included |
| Figure 2 | `analyses/cultural_pca.py`; `china/fit_china_mixed.m`; `analyses/resampling.py`: `snr_curve`; `analyses/map_statistics.py` | Cultural matrix, individual cortex and covariates, gradients, province residuals, spin indices |
| Figure 3 | `china/fit_china_relocation.m`; `analyses/map_statistics.py`: `commonality`; `analyses/group_summaries.py`: `shared_roi_bootstrap` | Birth/destination culture, mover status, ROI membership and participant-level ROI measures |
| Figure 4 | `china/fit_china_mediation.m`, `behavior_pcs.m`, `best_task_order.m`, `export_china_source.m` | Task metrics and directions, cortex, culture, covariates, gradient maps |
| Figure 5 | `international/prepare_study_cortex.py`, `prepare_study_units.py`, `fig5_chain.py`; `figures/render_fig5.m` | Raw country/state cortex, batch/demographic records, named macros, unit culture, aligned spatial targets and rotations; see [Figure 5 inputs](../international/FIGURE5.md) |
| Extended Data 1 | `analyses/psychometrics.py`; `analyses/cultural_pca.py`: `cronbach_alpha`; `analyses/group_summaries.py` | Item responses, grouping IDs, response-scale range; demographic and cortical measurements |
| Extended Data 2 | `workflows/study_inputs/generate_china_trajectory_inputs.m`, `generate_china_trajectory.py`, `assemble_center_geography.py`; `international/trajectory_models.py`; simulation and threshold entries | Raw cortex, center/province score orders, geography keys and international inputs; all-PC curves and geographic aggregates are generated |
| Extended Data 3 | `china/generate_gradient_inputs.m`; `workflows/run_specialized.py --analysis gradients` | Hemisphere ROI gradient summaries and surface eigenmodes; nine maps and correlation matrix |
| Extended Data 4 | `china/fit_china_mixed.m`; `analyses/harmonization.py`; `analyses/artifact_pca.py`; `analyses/plot_matrix_results.m` | Vertex/ROI cortical matrices, batch and biological designs, cortical mesh and spatial rotations; surface preprocessing is external |
| Extended Data 5 | `china/prepare_combat_inputs.m`; `workflows/study_inputs/run_combat_study_models.py`, `assemble_combat_panels.py`; shared mixed-model branch | Named raw cortical/covariate/culture fields, gradients and spins; protection settings saved separately |
| Extended Data 6 | `workflows/study_inputs/generate_culture_effect_rsa.m`, `render_culture_rsa.m` | Brain, raw residual, effect and gradient MAT inputs plus rotations; signed cultural RDM, observed RSA and direct-shuffle nulls are generated |
| Extended Data 7 | `workflows/study_inputs/run_area_volume.py`; `run_thickness_models.py`; `assemble_thickness_panels.py` | Named area/volume/thickness columns, PC3, gradient/eigenmode references and spins |
| Extended Data 8 | `china/run_model_branches.m`; `workflows/run_specialized.py --analysis rice_culture` / `rice_brain` | Rice share/suitability, cultural scores, thickness/area models, gradients and rotations |
| Extended Data 9 | `workflows/study_inputs/generate_china_resampling.m`, `export_china_resampling.py`; `international/resample_models.py` | Raw brain/residual inputs or the Fig.5 input folder; subject and whole-unit sampling designs, fit status, draws and quantiles are saved |
| Extended Data 10 | `china/run_china_deletion.m`, `render_china_deletion.m`; `international/run_loo.py` | Shared-chain inputs and outputs; harmonization/scaling and full-sample cultural residuals are fixed, artifacts and final models are refitted by fold |
| Extended Data 11 | `china/fit_china_relocation.m`; `china/run_sensitivity.m`; `analyses/group_summaries.py` | Relocation/cultural-distance variables, ROI groups and alternative adjustment columns |
| Extended Data 12 | `china/best_task_order.m`; `china/fit_china_mediation.m`; `analyses/map_statistics.py` | Seven task-effect vectors, behavioral PC effects, spatial reference data; one selected task order |
| Extended Data 13 | `workflows/study_inputs/run_thickness_models.py`, `assemble_dimension_correspondence.py`; `international/dimension_correspondence.py` | Raw cortex, covariates, six country/China dimensions, reference gradients and spins; LTO is long-term orientation |
| Extended Data 14 | `workflows/study_inputs/run_thickness_models.py --family hof6`, `generate_hofstede_webb.m`, `assemble_thickness_panels.py` | Alternative Hofstede cultural PCs, raw cortex, province influence, reference gradients and spins |
| Extended Data 15 | `analyses/cultural_pca.py`; `workflows/study_inputs/assemble_dimension_correspondence.py`; `international/dimension_correspondence.py` | Joint culture table and six computed dimension-effect maps; explicit PC orientations |
| Supplementary 1 | `analyses/run_word_pipeline.py`; `analyses/word_scores.py`; `analyses/cultural_pca.py`: `cronbach_alpha`; `analyses/map_statistics.py` | Raw daily keyword index workbook, province matching and external validity scores; cleaning and weighted PCA are included |
| Supplementary 2 | `analyses/traveling_calibration.py`: `fit_predict` | Repeated traveling-participant cortical measurements aligned across sites and ROIs; cross-participant train/test pairing must be supplied |
| Supplementary 3 | `simulation/generate_si3.m`, `export_si3_source.py`, `plot_si3.m` | Fully synthetic; complete reference realizations are included |
| Supplementary 4 | workflows/assemble.py; figures/layouts/reproduce_suppfig4.m | Selected participant relocation duration and birth/destination province labels |

`analyses/plot_matrix_results.m` exports scatter, matrix, mean/SEM and surface inputs to `.fig`. China and international examples have dedicated plotters. Original panel-specific map layouts and every empirical source-data exporter are not all implemented by this generic plotter.

Dedicated assembly contracts and .fig rendering: [workflows/README.md](../workflows/README.md). The entries above describe upstream requirements; the workflow README specifies which panel components are assembled automatically.

Reference meshes, geographic boundaries and nine dedicated original layouts are available via [figures/README.md](../figures/README.md).

Specialized assembly, original-layout rendering and fitted-result bridges: [workflows/SPECIALIZED_INPUTS.md](../workflows/SPECIALIZED_INPUTS.md).

Study-specific raw-input continuations, output files and plot entry points are listed in [study workflows](../workflows/study_inputs/README.md).


Additional leave-one-dimension-out sensitivity models use `china/run_dimension_omission.m`, `china/run_dimension_branches.m` and `international/dimension_maps.py`. Their omission plots are saved under `CulturalDimensionOmission`.
