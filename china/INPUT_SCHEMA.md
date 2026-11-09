# China model input

`T`: one row per participant, matched to the rows of `cortex` and all behavioral matrices. `cortex`: numeric participants × ROIs matrix of model-ready cortical measurements. `behavior`: cell array of seven participants × metrics matrices; task names are CVLT, GoNogo, NBack, Passball, RO, SDMT, UG.

Required table columns:

- The 28 numeric covariates returned by `china_covariates`.
- `province_birth`, `city_birth`: birthplace grouping IDs.
- `culture_pc1`, `culture_pc2`: province-level cultural scores mapped to participants.
- `relocator`: 0 for stayers, 1 for relocators.
- `destination_pc1`, `destination_pc2`: destination province cultural scores. Stayers use their birthplace scores.
- `province_destination`: destination province ID, required by `shuffle_china_culture`.

`prepare_china_table` constructs province-by-city IDs, standardizes covariates and birthplace culture, applies the same cultural scaling to destination culture, and sets cultural differences to zero for stayers. Use complete, finite input matrices with nonconstant covariate columns and a full-rank fixed-effects design. Sex uses a numeric binary encoding. Additional covariates should be scaled before model fitting.

The optional adjustments in Extended Data Table 1b can be supplied as columns and passed in `additional`: paddy-field percentage, rice suitability, province-mean genetic PC1 and PC2, human-transmissible disease rate, overall disease rate, cultivated-land proportion, clan density, and jinshi density. Their column names are caller-defined.

`behavior_pcs(behavior,directions)` accepts one vector per task, with +1 for higher-is-better metrics and −1 for lower-is-better metrics. Loadings, explained variance, retained component counts and directions are saved. Synthetic example metrics all use +1.
