clear; clc; close all;

if strcmp(getenv('V3_SMOKE'), '1')
    nSims_per_bin = 2;  out_name = 'sim_data_unified_v3_beta2x_matched_range6p0_smoke.mat';
    fprintf('SMOKE TEST MODE: nSims_per_bin = %d\n', nSims_per_bin);
else
    nSims_per_bin = 100; out_name = 'sim_data_unified_v3_beta2x_matched_range6p0.mat';
end
root=fileparts(fileparts(mfilename('fullpath')));
out_dir=fullfile(root,'outputs','simulation','ed2');
if strcmp(getenv('V3_SMOKE'),'1'),out_dir=fullfile(root,'outputs','simulation','ed2_smoke');end
if ~isfolder(out_dir),mkdir(out_dir);end

nSubjects = 5000; nROIs = 100; min_group_size = 30; noise_level = 0.5;
range_artif = 6.0; base_artif = 1.5; orth_thres = 0.005;      % artifact all-positive [1.50, 7.50], ~15-75x GT
mu_pos = 0.10; amp_pos = 0.03;           % all-positive GT: baseline + fluctuation (~0.1 effect size)

artif_rhos = [0, 0.1, 0.6, 0.8, 1.0];
target_rho_full_range = [-1, 1];
target_rho_null_range = [-0.1, 0.1];
nBins = length(artif_rhos);
nNull_per_bin = round(0.4 * nSims_per_bin);   % dense null-target sims
nTotal = nBins * nSims_per_bin;

rng(42,'twister');

out_cGT  = nan(1, nTotal);
out_sGTM = nan(1, nTotal);  out_eGTM = nan(1, nTotal);
out_sSAM = nan(1, nTotal);  out_eSAM = nan(1, nTotal);
out_sGTC = nan(1, nTotal);  out_eGTC = nan(1, nTotal);
out_sSAC = nan(1, nTotal);  out_eSAC = nan(1, nTotal);
out_maps_GT = nan(nTotal, nROIs);
out_maps_sM = nan(nTotal, nROIs);  out_maps_eM = nan(nTotal, nROIs);
out_maps_sC = nan(nTotal, nROIs);  out_maps_eC = nan(nTotal, nROIs);

is_null_sims = false(nBins, nSims_per_bin);
for ibin = 1:nBins
    is_null_sims(ibin, nSims_per_bin-nNull_per_bin+1:end) = true;
end

ideal_gradient = linspace(0.1, 0.3, nROIs)';

nProvinces_M = 30; nCenters_M = 29; max_pcs_M = 28;
nProvinces_C = 29; nCenters_C = 29; max_pcs_C = 27;

targG_M = zeros(nSubjects, 1);
counts = ones(nProvinces_M, 1) * min_group_size;
rem_subs = nSubjects - sum(counts);
props = rand(nProvinces_M, 1); props = props / sum(props);
counts = counts + floor(props * rem_subs);
counts(1) = counts(1) + (nSubjects - sum(counts));
idx = 1;
for i = 1:nProvinces_M
    targG_M(idx:idx+counts(i)-1) = i;
    idx = idx+counts(i);
end
targG_M = targG_M(randperm(nSubjects));

artifG_M = zeros(nSubjects, 1);
base_centers = 1:nCenters_M;
for i = 1:nProvinces_M
    subj_idx = find(targG_M == i);
    if i <= nCenters_M; main_center = i; else; main_center = randi(nCenters_M); end
    for j = 1:length(subj_idx)
        if rand() < 0.8
            artifG_M(subj_idx(j)) = main_center;
        else
            other_centers = base_centers(base_centers ~= main_center);
            artifG_M(subj_idx(j)) = other_centers(randi(length(other_centers)));
        end
    end
end

A_targ = zeros(nSubjects, 1);
counts = ones(nProvinces_C, 1) * floor(nSubjects / nProvinces_C);
rem_subs = nSubjects - sum(counts);
counts(1) = counts(1) + rem_subs;
idx = 1;
for i = 1:nProvinces_C
    A_targ(idx:idx+counts(i)-1) = i;
    idx = idx+counts(i);
end
A_targ = A_targ(randperm(nSubjects));
B_artf = A_targ;
artifG_C = B_artf; targG_C = A_targ;

shared_prov_vals = randn(nProvinces_M, 1);
real_prov_vals_M = shared_prov_vals(1:nProvinces_M);
real_prov_vals_C = shared_prov_vals(1:nProvinces_C);
vec_targ_real_M = zscore(real_prov_vals_M(targG_M));
vec_targ_real_C = zscore(real_prov_vals_C(targG_C));

is_valid_shared = false;
while ~is_valid_shared
    shared_center_noise = randn(nCenters_M, 1);
    vec_artif_M = zscore(shared_center_noise(artifG_M));
    vec_artif_C = zscore(shared_center_noise(artifG_C));
    is_valid_shared = abs(corr(vec_targ_real_M, vec_artif_M)) < orth_thres && ...
                      abs(corr(vec_targ_real_C, vec_artif_C)) < orth_thres;
end

shared_perm_scores = zscore(randn(nCenters_M, 3000));
perm_pool_M = shared_perm_scores(artifG_M, :);
perm_pool_C = shared_perm_scores(artifG_C, :);
clean_pool_real_M = perm_pool_M(:, abs(corr(perm_pool_M, vec_targ_real_M)) < orth_thres);
clean_pool_real_C = perm_pool_C(:, abs(corr(perm_pool_C, vec_targ_real_C)) < orth_thres);

Q_theory_M = orth(clean_pool_real_M - mean(clean_pool_real_M, 1));
rho_sq_theory_M = sum((Q_theory_M' * vec_targ_real_M).^2) / sum(vec_targ_real_M.^2);
k_theory_M = 1 - rho_sq_theory_M;
sigma_eff_theory_M = noise_level * sqrt((nSubjects - 1 - max_pcs_M) / (nSubjects - 1));
theory_file = fullfile(out_dir, 'theory_params_unified_v3_beta2x_matched_range6p0.mat');
save(theory_file, 'rho_sq_theory_M', 'k_theory_M', 'sigma_eff_theory_M', ...
    'noise_level', 'nSubjects', 'max_pcs_M');
fprintf('MIXED theory: joint R^2 = %.6f, rho = %.6f, k = %.6f, sigma_eff = %.6f\n', ...
    rho_sq_theory_M, sqrt(rho_sq_theory_M), k_theory_M, sigma_eff_theory_M);

if strcmp(getenv('V3_SETUP_ONLY'), '1')
    fprintf('SETUP-ONLY MODE: saved %s\n', theory_file);
    return;
end

fprintf('MIXED setup: corr(t,a) = %+.4f (target 0), pool = %dx%d\n', ...
    corr(vec_targ_real_M, vec_artif_M), size(clean_pool_real_M,1), size(clean_pool_real_M,2));
fprintf('COINC setup: corr(t,a) = %+.4f (target 0), pool = %dx%d\n', ...
    corr(vec_targ_real_C, vec_artif_C), size(clean_pool_real_C,1), size(clean_pool_real_C,2));

tasks_ibin = repmat((1:nBins)', 1, nSims_per_bin);   % column-major task table
tasks_isim = repmat(1:nSims_per_bin, nBins, 1);
tasks_ibin = tasks_ibin(:); tasks_isim = tasks_isim(:);

X_base_M = [ones(nSubjects,1), zscore(vec_targ_real_M)];
X_base_C = [ones(nSubjects,1), zscore(vec_targ_real_C)];

if ~strcmp(getenv('V3_SMOKE'), '1')
    nW = min(8, feature('numcores'));
    p = gcp('nocreate');
    if isempty(p)
        try
            parpool(nW);
        catch
            parpool('local', nW);
        end
    end
    fprintf('Parallel pool ready: %d workers\n, starting %d tasks\n', nW, nTotal);
end

t0 = tic;
max_workers=Inf;if strcmp(getenv('V3_SMOKE'),'1'),max_workers=0;end
parfor (itask = 1:nTotal,max_workers)
    ibin = tasks_ibin(itask); isim = tasks_isim(itask);
    rng(1000000 + ibin * 1000 + isim,'Threefry');

    artif_rho = artif_rhos(ibin);

    if is_null_sims(ibin, isim)
        target_rho = target_rho_null_range(1) + diff(target_rho_null_range) * rand(1);
    else
        target_rho = target_rho_full_range(1) + diff(target_rho_full_range) * rand(1);
    end

    spatial_noise = randn(nROIs, 1);
    spatial_noise = zscore(spatial_noise - ideal_gradient * (ideal_gradient \ spatial_noise));
    beta_GT_raw = mu_pos + amp_pos * (target_rho * zscore(ideal_gradient) + sqrt(1 - target_rho^2) * spatial_noise);
    beta_GT = max(beta_GT_raw, 0.02);                   % enforce strict positivity

    artif_noise = randn(nROIs, 1);
    artif_noise = zscore(artif_noise - ideal_gradient * (ideal_gradient \ artif_noise));
    raw_artif_shape = artif_rho * zscore(ideal_gradient) + sqrt(1 - artif_rho^2) * artif_noise;
    raw_artif_norm = (raw_artif_shape - min(raw_artif_shape)) / (max(raw_artif_shape) - min(raw_artif_shape));
    beta_Artif = base_artif + range_artif * raw_artif_norm;

    shared_noise_Y = randn(nSubjects, nROIs) * noise_level;

    Y_all_M = vec_targ_real_M * beta_GT' + vec_artif_M * beta_Artif' + shared_noise_Y;

    beta_est_start_M = zeros(nROIs, 1);
    beta_est_end_M = zeros(nROIs, 1);
    for r = 1:nROIs
        y_r = Y_all_M(:, r);
        b_raw = X_base_M \ zscore(y_r);
        beta_est_start_M(r) = b_raw(2);

        w_r = (clean_pool_real_M' * y_r) / (nSubjects - 1);
        [~, PCs_r, ~, ~, exp_r] = pca(clean_pool_real_M .* w_r', 'NumComponents', max_pcs_M);
        X_conf = [ones(nSubjects,1), PCs_r(:, 1:min(max_pcs_M, length(exp_r)))];
        res_y = y_r - X_conf * (X_conf \ y_r);
        b_final = X_base_M \ zscore(res_y);
        beta_est_end_M(r) = b_final(2);
    end

    Y_all_C = vec_targ_real_C * beta_GT' + vec_artif_C * beta_Artif' + shared_noise_Y;

    beta_est_start_C = zeros(nROIs, 1);
    beta_est_end_C = zeros(nROIs, 1);
    for r = 1:nROIs
        y_r = Y_all_C(:, r);
        b_raw = X_base_C \ zscore(y_r);
        beta_est_start_C(r) = b_raw(2);

        w_r = (clean_pool_real_C' * y_r) / (nSubjects - 1);
        [~, PCs_r, ~, ~, exp_r] = pca(clean_pool_real_C .* w_r', 'NumComponents', max_pcs_C);
        X_conf = [ones(nSubjects,1), PCs_r(:, 1:min(max_pcs_C, length(exp_r)))];
        res_y = y_r - X_conf * (X_conf \ y_r);
        b_final = X_base_C \ zscore(res_y);
        beta_est_end_C(r) = b_final(2);
    end

    out_cGT(itask)  = corr(beta_GT, ideal_gradient);
    out_sGTM(itask) = corr(beta_est_start_M, beta_GT);
    out_eGTM(itask) = corr(beta_est_end_M, beta_GT);
    out_sSAM(itask) = corr(beta_est_start_M, ideal_gradient);
    out_eSAM(itask) = corr(beta_est_end_M, ideal_gradient);
    out_sGTC(itask) = corr(beta_est_start_C, beta_GT);
    out_eGTC(itask) = corr(beta_est_end_C, beta_GT);
    out_sSAC(itask) = corr(beta_est_start_C, ideal_gradient);
    out_eSAC(itask) = corr(beta_est_end_C, ideal_gradient);
    out_maps_GT(itask, :) = beta_GT';
    out_maps_sM(itask, :) = beta_est_start_M';
    out_maps_eM(itask, :) = beta_est_end_M';
    out_maps_sC(itask, :) = beta_est_start_C';
    out_maps_eC(itask, :) = beta_est_end_C';
end
fprintf('All %d tasks done in %.1f min.\n', nTotal, toc(t0)/60);

corr_GT_SA = reshape(out_cGT, nBins, nSims_per_bin);
corr_est_start_GT_M = reshape(out_sGTM, nBins, nSims_per_bin);
corr_est_end_GT_M   = reshape(out_eGTM, nBins, nSims_per_bin);
corr_est_start_SA_M = reshape(out_sSAM, nBins, nSims_per_bin);
corr_est_end_SA_M   = reshape(out_eSAM, nBins, nSims_per_bin);
corr_est_start_GT_C = reshape(out_sGTC, nBins, nSims_per_bin);
corr_est_end_GT_C   = reshape(out_eGTC, nBins, nSims_per_bin);
corr_est_start_SA_C = reshape(out_sSAC, nBins, nSims_per_bin);
corr_est_end_SA_C   = reshape(out_eSAC, nBins, nSims_per_bin);

rho_ta_M = corr(vec_targ_real_M, vec_artif_M);
rho_ta_C = corr(vec_targ_real_C, vec_artif_C);
pool_size_M = size(clean_pool_real_M);
pool_size_C = size(clean_pool_real_C);

maps_GT = out_maps_GT;
maps_sM = out_maps_sM;  maps_eM = out_maps_eM;
maps_sC = out_maps_sC;  maps_eC = out_maps_eC;

save(fullfile(out_dir, out_name), 'artif_rhos', 'nSims_per_bin', 'nNull_per_bin', ...
    'corr_GT_SA', 'is_null_sims', 'target_rho_full_range', 'target_rho_null_range', ...
    'mu_pos', 'amp_pos', 'nProvinces_M', 'nCenters_M', 'max_pcs_M', ...
    'nProvinces_C', 'nCenters_C', 'max_pcs_C', ...
    'rho_ta_M', 'rho_ta_C', 'pool_size_M', 'pool_size_C', ...
    'corr_est_start_GT_M', 'corr_est_end_GT_M', 'corr_est_start_SA_M', 'corr_est_end_SA_M', ...
    'corr_est_start_GT_C', 'corr_est_end_GT_C', 'corr_est_start_SA_C', 'corr_est_end_SA_C', ...
    'maps_GT', 'maps_sM', 'maps_eM', 'maps_sC', 'maps_eC');
disp(['Simulation completed and data saved to ' out_name]);
