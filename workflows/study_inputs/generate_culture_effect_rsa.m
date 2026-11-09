function generate_culture_effect_rsa(config_file,output_dir)
% Generate province RDMs, cultural RSA effects and fixed-seed permutation results.
cfg=jsondecode(fileread(config_file));if ~isfolder(output_dir),mkdir(output_dir);end
rot=load(cfg.spins_file,'perm_id');
if ~isfield(cfg,'seed'),cfg.seed=20261029;end
if ~isfield(cfg,'surface_statistic'),cfg.surface_statistic='pearson';end
if ~isfield(cfg,'rsa_statistic'),cfg.rsa_statistic='spearman';end
assert(ismember(cfg.rsa_statistic,{'spearman','pearson'}));
assert(ismember(cfg.surface_statistic,{'spearman','pearson'}));rng(cfg.seed,'twister');

fprintf('Loading data and initializing...\n');

rot_times = min(1000,size(rot.perm_id,2));
load(cfg.dataset_file);
load(cfg.brain_file,'brain_data');

addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
controlvars=china_covariates;
all_vars=fieldnames(brain_data);
brain_var=all_vars(contains(all_vars,'thickness') & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}));
S=load(cfg.raw_file,'raw_res','valididx');
assert(isfield(S,'valididx'),'Residual input must include its original valididx row mask');
valididx=S.valididx(:);assert(numel(valididx)==numel(brain_data.age));
T=study_china_table(brain_data,find(valididx));
brainM_valid=S.raw_res{1,29};assert(size(brainM_valid,1)==height(T));
final_valid_idx=valididx;
brain_data.pc7_2=-brain_data.pc7_2;
pc1_valid=brain_data.pc7_1(valididx);pc2_valid=brain_data.pc7_2(valididx);
brainM_clean=study_china_adjust(T,brainM_valid,[pc1_valid pc2_valid]);
n_subjects=height(T);n_rois=size(brainM_clean,2);
valididx_rsa = final_valid_idx;
prov_vec = brain_data.province_birth(valididx_rsa);
unique_prov = unique(brain_data.province_birth(valididx_rsa));
cohend = @(v1,v2) (mean(v1) - mean(v2));

RDM_brain = nan(length(unique_prov),length(unique_prov),n_rois);
for roi = 1:n_rois
    for cy1 = 1:length(unique_prov)
        cy1idx = prov_vec == unique_prov(cy1);
        for cy2 = 1:length(unique_prov)
            cy2idx = prov_vec == unique_prov(cy2);
            RDM_brain(cy1,cy2,roi) = cohend(brainM_clean(cy1idx,roi),brainM_clean(cy2idx,roi));
        end
    end
end

Culnames = {'pc7_1','pc7_2'};
RDM_culture = nan(length(unique_prov),length(unique_prov),length(Culnames));
for ccc = 1:length(Culnames)
    cul_vec = brain_data.(Culnames{ccc})(valididx_rsa);
    for cy1 = 1:length(unique_prov)
        cy1idx = prov_vec == unique_prov(cy1);
        for cy2 = 1:length(unique_prov)
            cy2idx = prov_vec == unique_prov(cy2);
            RDM_culture(cy1,cy2,ccc) = mean(cul_vec(cy1idx)) - mean(cul_vec(cy2idx));
        end
    end
end

n_prov = length(unique_prov);
r_matrix = zeros(length(Culnames), n_rois);
for k = 1:length(Culnames)
    cul_diff_k = RDM_culture(:, :, k);
    cul_diff_vec = cul_diff_k(triu(true(n_prov), 1));
    for j = 1:n_rois
        cohend_diff_j = RDM_brain(:, :, j);
        cohend_diff_vec = cohend_diff_j(triu(true(n_prov), 1));
        [r, p] = corr(cul_diff_vec, cohend_diff_vec,'Type', cfg.rsa_statistic);
        r_matrix(k, j) = r;
    end
end

fprintf('Starting province permutation test...\n');
n_permutations = 10000;
if isfield(cfg,'draws'),n_permutations=cfg.draws;end
null_distribution = zeros(2, n_permutations, n_rois);accepted_indices=zeros(2,n_permutations,n_prov);spin_rng_states=cell(2,1);

mask = triu(true(n_prov), 1);
upper_tri_indices = find(mask);
n_pairs = length(upper_tri_indices);
brain_vectors = zeros(n_pairs, n_rois);
for roi = 1:n_rois
    temp_matrix = RDM_brain(:, :, roi);
    brain_vectors(:, roi) = temp_matrix(upper_tri_indices);
end

for ccc = 1:2
    Cul_matrix = squeeze(RDM_culture(:,:,ccc));
    cul_vector_observed = Cul_matrix(upper_tri_indices);

    for p = 1:n_permutations
        perm=randperm(n_prov);
        permuted_Cul=Cul_matrix(perm,perm);
        vec_cul_perm=permuted_Cul(upper_tri_indices);
        accepted_indices(ccc,p,:)=perm;null_distribution(ccc, p, :) = corr(brain_vectors, vec_cul_perm,'Type',cfg.rsa_statistic);
    end
end

save(fullfile(output_dir,'RSA_checkpoint.mat'),'r_matrix','null_distribution','accepted_indices','RDM_brain','RDM_culture','-v7.3');
fprintf('Processing display components...\n');
load(cfg.effect_file,'AllRes');
CulturePC1 = AllRes.Expcs.Coeff(:,1);
CulturePC2 = AllRes.Expcs.Coeff(:,2);

ED5_bc_Data = struct();
ED5_bc_Data.CulturePC1 = CulturePC1;
ED5_bc_Data.CulturePC2_raw = CulturePC2; % Will invert manually during plotting or here

for ccc = 1:2
    rsa_obs_r = r_matrix(ccc,:)';
    null_dist_ccc = squeeze(null_distribution(ccc, :, :));
    surface_observed=rsa_obs_r;
    z_val = (surface_observed - mean(null_dist_ccc, 1)') ./ std(null_dist_ccc, 0, 1)';
    climsss = ceil(max(abs(z_val)) * 1.15) * [-1 1];

    if ccc == 1
        brainmap_tgt = -mode2;
        hist_ref_name = 'Eigenmode2';
        scatter_y = zscore(CulturePC1);
    else
        brainmap_tgt = SAaxis;
        hist_ref_name = 'S-A Axis';
        scatter_y = -zscore(CulturePC2);
    end

    tpR = corr(null_dist_ccc', brainmap_tgt(:), 'type',cfg.surface_statistic);
    real_value = corr(rsa_obs_r, brainmap_tgt(:),'type',cfg.surface_statistic);

    scatter_x = zscore(rsa_obs_r);
    spin_rng_states{ccc}=rng; [p_values_spin,spin_null] = perm_sphere_p(scatter_x,scatter_y,rot.perm_id(:,1:rot_times),'pearson');

    if ccc == 1
        ED5_bc_Data.PC1.brainmap_z = z_val;
        ED5_bc_Data.PC1.climsss = climsss;
        ED5_bc_Data.PC1.hist_real_val = real_value;
        ED5_bc_Data.PC1.hist_null_dist = tpR;
        ED5_bc_Data.PC1.hist_ref_name = hist_ref_name;
        ED5_bc_Data.PC1.scatter_x = scatter_x;
        ED5_bc_Data.PC1.scatter_y = scatter_y;
        ED5_bc_Data.PC1.spin_p = p_values_spin;
        ED5_bc_Data.PC1.spin_null = spin_null;
    else
        ED5_bc_Data.PC2.brainmap_z = z_val;
        ED5_bc_Data.PC2.climsss = climsss;
        ED5_bc_Data.PC2.hist_real_val = real_value;
        ED5_bc_Data.PC2.hist_null_dist = tpR;
        ED5_bc_Data.PC2.hist_ref_name = hist_ref_name;
        ED5_bc_Data.PC2.scatter_x = scatter_x;
        ED5_bc_Data.PC2.scatter_y = scatter_y;
        ED5_bc_Data.PC2.spin_p = p_values_spin;
        ED5_bc_Data.PC2.spin_null = spin_null;
    end
end

save(fullfile(output_dir,'RSA_results.mat'),'ED5_bc_Data','r_matrix','null_distribution','accepted_indices','spin_rng_states','RDM_brain','RDM_culture','-v7.3');
save(fullfile(output_dir,'RSA_source_results.mat'),'ED5_bc_Data','cfg','-v7');
fprintf('Source Data processing completed and saved successfully!\n');
save(fullfile(output_dir,'generation_inputs.mat'),'cfg','valididx','brain_var','unique_prov');
end

