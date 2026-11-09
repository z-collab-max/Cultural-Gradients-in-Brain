# China individual-model inputs

The individual cortical models use `china/china_covariates.m` and `china/china_formula.m`: 28 fixed covariates, a predictor, an intercept, and random intercepts for birthplace province and city within province. Fits use MATLAB `fitlme` with REML. Each permutation refits the same model, including its variance components. Scanner adjustment, ComBat and artifact construction remain separate preprocessing steps. Spatial regressions between ROI maps and psychometric calculations retain their own definitions.

Python model adapters require MATLAB and the Statistics and Machine Learning Toolbox. Put `matlab` on PATH or set `MATLAB_EXECUTABLE` to your local MATLAB executable. Run commands from the repository root. The Python adapter calls the same MATLAB formula as the native workflows; it does not substitute an OLS fit. A full 10,000-draw run over hundreds of ROIs can take substantial time; use a small `--draws` value to check input alignment first.

Supply these numeric vectors in `brain_data`, in the same participant order as the cortical measurements:

| Fields | Meaning |
|---|---|
| `age`, `sex`, `education` | Age, numeric sex coding, education |
| `genetic_pc1` … `genetic_pc10` | Ten genetic principal components |
| `father_education`, `mother_education` | Parental education |
| `intracranial_volume`, `total_surface_area` | Global cortical covariates |
| `absolute_latitude`, `terrain_ruggedness` | Absolute latitude and terrain ruggedness |
| `historical_humidity`, `historical_temperature`, `historical_sunshine`, `historical_precipitation` | Historical climate |
| `population_density`, `green_space` | Population density and green space |
| `natural_disaster_frequency_per_capita` | Per-capita natural-disaster frequency |
| `hdi_education`, `hdi_income` | HDI education and income components |

Also supply `province_birth` and `province_city` as finite numeric group codes. The city code identifies birthplace city; the model combines it with province so repeated local city codes in different provinces are distinct. Do not use recruitment site as birthplace city. Canonical field names are required: callers must explicitly map their local fields and units before running. Missing fields, constant covariates and rank-deficient fixed designs are errors; there is no reduced-covariate fallback. Continuous columns and predictors are sample-SD standardized before fitting.

Retain the outcome, culture and scanner fields specified by each workflow. Province permutation indices must follow sorted `province_birth` codes and contain each province once per draw. Cultural predictors must be constant within province. Supplied culture tables must use the same province order. No draws are screened against the original culture scores.

`china/prepare_combat_inputs.m` exports the 28 fields and both group codes under `model_design` in the ComBat MAT file. `run_combat_study_models.py` requires that structure. Old intermediate files without `model_design` must be regenerated. ComBat protection settings are independent preprocessing choices; observed and null final fits use the same protected matrix and the same mixed-model design.

Residual MAT files for RSA and resampling require `raw_res` and the original full-length logical `valididx` mask, so row identity is preserved. Resampling records unsuccessful fits in `fit_failures` and leaves their maps as NaN; downstream summaries must use valid-fit counts.

The public examples use supplied or synthetic inputs. These model examples are not a promise of numerical identity with published empirical results.
