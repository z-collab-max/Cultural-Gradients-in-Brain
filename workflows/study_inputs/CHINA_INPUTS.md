# China individual-model inputs

The individual cortical models use `china/china_covariates.m` and `china/china_formula.m`: 17 individual covariates plus all nonzero city/province macro PCs, a predictor, an intercept, and random intercepts for birthplace province and city within province. Fits use MATLAB `fitlme` with REML. Each permutation refits the same model, including its variance components. Scanner adjustment, ComBat and artifact construction remain separate preprocessing steps. Spatial regressions between ROI maps and psychometric calculations retain their own definitions.

Python model adapters require MATLAB and the Statistics and Machine Learning Toolbox. Put `matlab` on PATH or set `MATLAB_EXECUTABLE` to your local MATLAB executable. Run commands from the repository root. The Python adapter calls the same MATLAB formula as the native workflows; it does not substitute an OLS fit. A full 10,000-draw run over hundreds of ROIs can take substantial time; use a small `--draws` value to check input alignment first.

First run the separate [macro PCA preprocessing](../../covariates/README.md). Use `china/prepare_brain_macro_pcs.m` for a `brain_data` MAT file. Supply these numeric vectors in `brain_data`, in the same participant order as the cortical measurements:

| Fields | Meaning |
|---|---|
| `age`, `sex`, `education` | Age, numeric sex coding, education |
| `genetic_pc1` … `genetic_pc10` | Ten genetic principal components |
| `father_education`, `mother_education` | Parental education |
| `intracranial_volume`, `total_surface_area` | Global cortical covariates |
| `city_macro_pc1` … | All nonzero PCs of city-level macro covariates |
| `province_macro_pc1` … | All nonzero PCs of province-level macro covariates |

Also supply `province_birth` and `province_city` as finite numeric group codes. The city code identifies birthplace city; the model combines it with province so repeated local city codes in different provinces are distinct. Do not use recruitment site as birthplace city. Canonical field names are required: callers must explicitly map their local fields and units before running. Missing fields, constant covariates and rank-deficient fixed designs are errors; there is no reduced-covariate fallback. Continuous columns and predictors are sample-SD standardized before fitting.

Retain the outcome, culture and scanner fields specified by each workflow. Province permutation indices must follow sorted `province_birth` codes and contain each province once per draw. Cultural predictors must be constant within province. Supplied culture tables must use the same province order. No draws are screened against the original culture scores.

`china/prepare_combat_inputs.m` exports the individual and macro-PC fields and both group codes under `model_design` in the ComBat MAT file. `run_combat_study_models.py` requires that structure. Old intermediate files without `model_design` must be regenerated. ComBat protection settings are independent preprocessing choices; observed and null final fits use the same protected matrix and the same mixed-model design.

Residual MAT files for RSA and resampling require `raw_res` and the original full-length logical `valididx` mask, so row identity is preserved. Resampling records unsuccessful fits in `fit_failures` and leaves their maps as NaN; downstream summaries must use valid-fit counts.

The examples use caller-supplied or synthetic inputs.
