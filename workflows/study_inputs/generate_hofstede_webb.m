function generate_hofstede_webb(brain_file,gradient_file,output_file,n_boot)
if nargin<4,n_boot=10000;end
assert(n_boot>0 && mod(n_boot,2)==0);
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
out=fileparts(output_file);if ~isempty(out)&&~isfolder(out),mkdir(out);end
rng(42, 'twister');

fprintf('========================================================================\n');
fprintf('  Wild Webb Bootstrap 95%% CI Calculation (Culture PCs & Gradients)\n');
fprintf('========================================================================\n\n');
fprintf('1. Loading 6 canonical cortical gradients...\n');
load(gradient_file, 'mode2', 'SAaxis', 'reardon_dev', 'xu_evo', 'intervar', 'cogPC1');
grad_names = {'eigenmode2', 'SAaxis', 'devexp', 'evoexp', 'intersubjvar', 'cogpc1'};
grad_matrix = [-mode2(:), SAaxis(:), reardon_dev(:), xu_evo(:), intervar(:), cogPC1(:)]; % [308 x 6]
Z_grad = zscore(grad_matrix);
n_rois = size(grad_matrix,1);
n_grads = 6;
Ridge_Projector = (Z_grad' * Z_grad) \ Z_grad'; % [6 x 308]
fprintf('2. Loading subject data and setting up China models...\n');

load(brain_file,'brain_data');


controlvars=china_covariates;
Culnames_birth_pcs = {'ValidP_hofPC6_1_birth','ValidP_hofPC6_2_birth','ValidP_hofPC6_3_birth'};
Culnames_recruit_pcs={};
all_vars = fieldnames(brain_data);
brainidx = (contains(all_vars,'thickness') & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}));
brain_var = all_vars(brainidx);
controlM = nan(length(brain_data.age),length(controlvars));
for cvv = 1:length(controlvars), controlM(:,cvv) = brain_data.(controlvars{cvv})(:); end
brainM = nan(length(brain_data.age),length(brain_var));
for bvv = 1:length(brain_var), brainM(:,bvv) = brain_data.(brain_var{bvv}); end
predM_birth = nan(length(brain_data.age),length(Culnames_birth_pcs));
for hvv = 1:length(Culnames_birth_pcs), predM_birth(:,hvv) = brain_data.(Culnames_birth_pcs{hvv}); end
predM_recruit = nan(length(brain_data.age),length(Culnames_recruit_pcs));
for hvv = 1:length(Culnames_recruit_pcs), predM_recruit(:,hvv) = brain_data.(Culnames_recruit_pcs{hvv}); end

valididx = ~isnan(sum(controlM,2)) & ~isnan(sum(predM_birth,2)) & ~isnan(sum(predM_recruit,2)) & ~isnan(sum(brainM,2)) & ~isnan(brain_data.centerID_rerange);
n_subs = sum(valididx);

center_code = brain_data.centerID_rerange(valididx);
centeriddd = dummyvar(categorical(center_code));
Y_expcs = nan(n_subs, n_rois);
for bvv = 1:n_rois
    [~,~,stats] = glmfit(centeriddd(:,1:end-1), brainM(valididx,bvv));
    Y_expcs(:,bvv) = zscore(stats.resid);
end

T=study_china_table(brain_data,find(valididx));
Province_cat_cn = brain_data.province_birth(valididx);
[prov_list_cn, ~, prov_idx_cn] = unique(Province_cat_cn);
n_prov = length(prov_list_cn);
fprintf('  Valid subjects: %d, Birth provinces: %d\n', n_subs, n_prov);
target_names = {
    'ValidP_hofPC6_1_birth', ...
    'ValidP_hofPC6_2_birth', ...
    'Hofstede_indvres_control', ...
    'Hofstede_powerdis_control'
};
target_display = {
    'Culture PC1 (6-dim Hofstede)', ...
    'Culture PC2 (6-dim Hofstede)', ...
    'IVR (Col 1 Hofstede indvres)', ...
    'PDI (Col 2 Hofstede powerdis)'
};
n_targets = length(target_names);
Obs = struct();
Models = struct();

for t = 1:n_targets
    t_name = target_names{t};
    raw_val = brain_data.(t_name)(valididx);
    X_target = zscore(raw_val);
    fitted=study_china_influence(T,Y_expcs,X_target);b_val=fitted.point;
    r_em2_raw = corr(b_val, -mode2(:));
    r_sa_raw  = corr(b_val, SAaxis(:));
    if (r_em2_raw < 0 || r_sa_raw < 0)
        sidir = -1;
    else
        sidir = +1;
    end
    b_aligned = sidir * b_val;
    r_em2 = corr(b_aligned, -mode2(:));
    r_sa  = corr(b_aligned, SAaxis(:));
    Z_y = zscore(b_aligned);
    multi_betas = Ridge_Projector * Z_y; % [6 x 1]
    
    Obs(t).name = t_name;
    Obs(t).display = target_display{t};
    Obs(t).sidir = sidir;
    Obs(t).b_aligned = b_aligned;
    Obs(t).r_em2 = r_em2;
    Obs(t).r_sa  = r_sa;
    Obs(t).multi_betas = multi_betas;
    
    Models(t).influence=fitted.influence;
    Models(t).b_aligned = b_aligned;
    Models(t).sidir = sidir;
end

fprintf('\n3. Running Level A: Province-Clustered Wild Webb Bootstrap (%d draws)...\n', n_boot);
webb = [-sqrt(1.5), -1, -sqrt(.5), sqrt(.5), 1, sqrt(1.5)];
half = n_boot / 2;
base_w = repmat(webb(:), ceil(half / 6), 1);
base_w = base_w(1:half);

mult_half_cn = zeros(half, n_prov);
for g = 1:n_prov
    mult_half_cn(:, g) = base_w(randperm(half));
end
multG_all_cn = [mult_half_cn; -mult_half_cn]; % [n_boot x n_prov] Balanced-antithetic

LevelA = struct();

for t = 1:n_targets
    fprintf('  Processing Target %d/%d: %s ...\n', t, n_targets, target_names{t});
    influence=Models(t).influence;
    Delta=multG_all_cn*influence;
    influence_all(t,:,:)=influence;
    Boot_Maps = Models(t).b_aligned' + Models(t).sidir * Delta;
    r_em2_boot = zeros(n_boot, 1);
    r_sa_boot  = zeros(n_boot, 1);
    
    target_em2 = -mode2(:);
    target_sa  = SAaxis(:);
    z_em2 = zscore(target_em2);
    z_sa  = zscore(target_sa);
    Z_Boot = zscore(Boot_Maps, 0, 2); % [n_boot x 308] standardized across ROIs
    r_em2_boot = (Z_Boot * z_em2) / (n_rois - 1);
    r_sa_boot  = (Z_Boot * z_sa) / (n_rois - 1);
    Multi_Betas_Boot = (Ridge_Projector * Z_Boot')'; % [n_boot x 6]
    
    DrawsA(t).r_em2=r_em2_boot;DrawsA(t).r_sa=r_sa_boot;DrawsA(t).multi_beta=Multi_Betas_Boot;
    LevelA(t).name = target_names{t};
    LevelA(t).display = target_display{t};
    LevelA(t).obs_r_em2 = Obs(t).r_em2;
    LevelA(t).obs_r_sa  = Obs(t).r_sa;
    LevelA(t).obs_multi = Obs(t).multi_betas;
    
    LevelA(t).r_em2_ci = prctile(r_em2_boot, [2.5, 97.5]);
    LevelA(t).r_em2_se = std(r_em2_boot);
    LevelA(t).r_em2_mean = mean(r_em2_boot);
    LevelA(t).r_em2_p = 2 * min(mean(r_em2_boot <= 0), mean(r_em2_boot >= 0));
    
    LevelA(t).r_sa_ci = prctile(r_sa_boot, [2.5, 97.5]);
    LevelA(t).r_sa_se = std(r_sa_boot);
    LevelA(t).r_sa_mean = mean(r_sa_boot);
    LevelA(t).r_sa_p = 2 * min(mean(r_sa_boot <= 0), mean(r_sa_boot >= 0));
    
    LevelA(t).multi_ci = prctile(Multi_Betas_Boot, [2.5, 97.5], 1); % [2 x 6]
    LevelA(t).multi_se = std(Multi_Betas_Boot, 0, 1)';               % [6 x 1]
    LevelA(t).multi_mean = mean(Multi_Betas_Boot, 1)';               % [6 x 1]
    LevelA(t).multi_p = zeros(n_grads, 1);
    for g = 1:n_grads
        LevelA(t).multi_p(g) = 2 * min(mean(Multi_Betas_Boot(:, g) <= 0), mean(Multi_Betas_Boot(:, g) >= 0));
    end
end
fprintf('\n4. Running Level B: Cortical ROI-Level Wild Webb Bootstrap (%d draws)...\n', n_boot);

mult_half_roi = zeros(half, n_rois);
for r = 1:n_rois
    mult_half_roi(:, r) = base_w(randperm(half));
end
mult_all_roi = [mult_half_roi; -mult_half_roi]; % [n_boot x 308]

LevelB = struct();

for t = 1:n_targets
    fprintf('  Processing Target %d/%d: %s ...\n', t, n_targets, target_names{t});
    
    Z_y = zscore(Obs(t).b_aligned); % [308 x 1]
    X_em2 = [ones(n_rois, 1), zscore(target_em2)];
    Ainv_em2 = pinv(X_em2' * X_em2);
    B_em2 = Ainv_em2 * (X_em2' * Z_y);
    resid_em2 = Z_y - X_em2 * B_em2;
    lev_em2 = sum((X_em2 * Ainv_em2) .* X_em2, 2);
    resid_em2_hc2 = resid_em2 ./ sqrt(max(1 - lev_em2, 1e-7));
    slope_em2 = Ainv_em2(2, :) * X_em2'; % [1 x 308]
    
    delta_em2 = (mult_all_roi .* slope_em2) * resid_em2_hc2;
    r_em2_boot = B_em2(2) + delta_em2;
    X_sa = [ones(n_rois, 1), zscore(target_sa)];
    Ainv_sa = pinv(X_sa' * X_sa);
    B_sa = Ainv_sa * (X_sa' * Z_y);
    resid_sa = Z_y - X_sa * B_sa;
    lev_sa = sum((X_sa * Ainv_sa) .* X_sa, 2);
    resid_sa_hc2 = resid_sa ./ sqrt(max(1 - lev_sa, 1e-7));
    slope_sa = Ainv_sa(2, :) * X_sa'; % [1 x 308]
    
    delta_sa = (mult_all_roi .* slope_sa) * resid_sa_hc2;
    r_sa_boot = B_sa(2) + delta_sa;
    Ainv_multi = pinv(Z_grad' * Z_grad);
    B_multi = Ainv_multi * (Z_grad' * Z_y); % [6 x 1]
    resid_multi = Z_y - Z_grad * B_multi;
    lev_multi = sum((Z_grad * Ainv_multi) .* Z_grad, 2);
    resid_multi_hc2 = resid_multi ./ sqrt(max(1 - lev_multi, 1e-7));
    Multi_Betas_Boot = zeros(n_boot, n_grads);
    for g = 1:n_grads
        slope_row_g = Ainv_multi(g, :) * Z_grad'; % [1 x 308]
        delta_g = (mult_all_roi .* slope_row_g) * resid_multi_hc2;
        Multi_Betas_Boot(:, g) = B_multi(g) + delta_g;
    end
    
    DrawsB(t).r_em2=r_em2_boot;DrawsB(t).r_sa=r_sa_boot;DrawsB(t).multi_beta=Multi_Betas_Boot;
    LevelB(t).name = target_names{t};
    LevelB(t).display = target_display{t};
    LevelB(t).obs_r_em2 = Obs(t).r_em2;
    LevelB(t).obs_r_sa  = Obs(t).r_sa;
    LevelB(t).obs_multi = Obs(t).multi_betas;
    
    LevelB(t).r_em2_ci = prctile(r_em2_boot, [2.5, 97.5]);
    LevelB(t).r_em2_se = std(r_em2_boot);
    LevelB(t).r_em2_mean = mean(r_em2_boot);
    LevelB(t).r_em2_p = 2 * min(mean(r_em2_boot <= 0), mean(r_em2_boot >= 0));
    
    LevelB(t).r_sa_ci = prctile(r_sa_boot, [2.5, 97.5]);
    LevelB(t).r_sa_se = std(r_sa_boot);
    LevelB(t).r_sa_mean = mean(r_sa_boot);
    LevelB(t).r_sa_p = 2 * min(mean(r_sa_boot <= 0), mean(r_sa_boot >= 0));
    
    LevelB(t).multi_ci = prctile(Multi_Betas_Boot, [2.5, 97.5], 1); % [2 x 6]
    LevelB(t).multi_se = std(Multi_Betas_Boot, 0, 1)';               % [6 x 1]
    LevelB(t).multi_mean = mean(Multi_Betas_Boot, 1)';               % [6 x 1]
    LevelB(t).multi_p = zeros(n_grads, 1);
    for g = 1:n_grads
        LevelB(t).multi_p(g) = 2 * min(mean(Multi_Betas_Boot(:, g) <= 0), mean(Multi_Betas_Boot(:, g) >= 0));
    end
end
save_path = output_file;
save(save_path, 'Obs', 'LevelA', 'LevelB', 'grad_names','DrawsA','DrawsB','multG_all_cn','mult_all_roi','prov_list_cn','influence_all','-v7');
fprintf('\nSaved results to: %s\n', save_path);


end
