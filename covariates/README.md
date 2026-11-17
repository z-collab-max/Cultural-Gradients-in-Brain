# Geographic macro-covariate PCA

Run this preprocessing step before cortical modeling. Cultural scores, individual covariates, brain measurements and outcomes are not included in this PCA.

| Cohort | PCA rows | Model columns | Full-rank count |
|---|---|---|---|
| China city | One row per birthplace province + city | `city_macro_pc1` … | 9 |
| China province | One row per birthplace province | `province_macro_pc1` … | 2 |
| USA | One row per state | `state_macro_pc1` … | 8 |
| Nations | One row per country | `nation_macro_pc1` … | 12 |

Each geographic unit receives equal weight. Raw columns are centred and sample-SD standardized across units. SVD retains **all nonzero components** using `max(matrix.shape) * machine_epsilon * largest_singular_value` as the numerical rank tolerance. Constant columns contribute no components. The largest absolute loading in each component is positive to give reproducible signs. Input rows must be complete, and macro values must be constant within each declared unit. At least one nonzero component per required level is needed by the model adapters.

The raw-to-PC field mapping is defined in `input_fields.py`; these fields belong in preprocessing inputs, not final model formulas. Means, scales, loadings, explained variance, rank, geographic keys and PC names are saved alongside prepared tables. Retaining all nonzero components preserves the span of the original nuisance variables together with the intercept. It does not guarantee identical results if samples, covariates, random effects or other preprocessing steps change.

## CSV inputs

From the repository root:

```sh
python covariates/prepare_macro_pcs.py --cohort China --input inputs/china/covariates_raw.csv --output inputs/china/covariates.csv
python covariates/prepare_macro_pcs.py --cohort USA --input inputs/international/USA/units_raw.csv --output inputs/international/USA/units.csv
python covariates/prepare_macro_pcs.py --cohort Nations --input inputs/international/Nations/units_raw.csv --output inputs/international/Nations/units.csv
```

China raw CSV contains `participant_id`, `province_birth`, `city_birth`, the 17 individual covariates and the 11 raw macro fields. The output preserves row order and individual fields, replaces the macros with PC columns and retains the grouping keys. International raw tables have unique `unit` rows, raw macro fields, cultural scores and optional geographic display labels. PC tables keep these identifiers and cultural scores. Each output CSV has a companion `.pca.json` file. Downstream models read the PC columns in numerical order and accommodate a reduced nonzero rank.

## MATLAB inputs

```matlab
addpath('china');
[T, macro_pca] = prepare_macro_pcs(T_raw, fullfile('inputs','china','prepared_table.mat'));
% Or prepare an existing brain_data structure for study adapters:
prepare_brain_macro_pcs(fullfile('inputs','china','brain_raw.mat'), ...
                       fullfile('inputs','china','brain_pc.mat'));
```

`T_raw` uses `province_birth` and `city_birth`; `brain_data` uses `province_birth` and `province_city`. Both create the same city-within-province keys. All original cortical, cultural and participant-ID fields in `brain_data` are retained in the local prepared file. Use that file as the study adapter's brain input. Keep individual input and intermediate files outside the distributed code.

The synthetic China and international input generators call these separate PCA functions automatically and save their PCA metadata. With empirical inputs, fit the PCA once on the declared set of geographic units. Apply saved means, scales and loadings when mapping to participants or held-out data; do not refit using outcomes, cultural associations or individual sampling frequencies. `macro_pca.transform` applies saved Python metadata. MATLAB scores are `((raw - metadata.mean) ./ metadata.scale) * metadata.loadings`.

Optional macro nuisance blocks for sensitivity analyses can be passed to Python `prepare(table, blocks)` with explicit geographic keys, raw fields and a distinct PC prefix, then supplied as `additional` model terms. Keep indicators whose cortical associations are being tested as predictors separate from nuisance-PC construction.
