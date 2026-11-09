clear; clc; close all;

if strcmp(getenv('V3_SMOKE'), '1')
    nSims_per_bin = 2;  out_name = 'sim_data_unified_v3_beta2x_matched_range6p0_smoke.mat';
    fprintf('SMOKE TEST MODE: nSims_per_bin = %d\n', nSims_per_bin);
else
    nSims_per_bin = 100; out_name = 'sim_data_unified_v3_beta2x_matched_range6p0.mat';
end
package=fileparts(fileparts(mfilename('fullpath')));
out_dir=fullfile(package,'outputs','simulation','si3');
if strcmp(getenv('V3_SMOKE'),'1')
    out_dir=fullfile(package,'outputs','simulation','si3_smoke');
end
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

rng(20260917,'twister');
shared_perm_scores = zscore(randn(nCenters_M, 3000));
perm_pool_M = shared_perm_scores(artifG_M, :);
perm_pool_C = shared_perm_scores(artifG_C, :);
clean_pool_real_M = perm_pool_M(:, abs(corr(perm_pool_M, vec_targ_real_M)) < orth_thres);
clean_pool_real_C = perm_pool_C(:, abs(corr(perm_pool_C, vec_targ_real_C)) < orth_thres);


source_counts=[1 3 5]; max_sources=max(source_counts);
constraints=accumarray(artifG_C,vec_targ_real_C,[nCenters_C,1]);
Q=orth(constraints);
rng(20260914,'twister');
center_sources_raw=randn(nCenters_M,max_sources);
center_sources_M=center_sources_raw;
center_sources_C=center_sources_raw-Q*(Q'*center_sources_raw);
scores_M=zscore(center_sources_M(artifG_M,:));
scores_C=zscore(center_sources_C(artifG_C,:));
source_culture_r_M=corr(vec_targ_real_M,scores_M);
source_culture_r_C=corr(vec_targ_real_C,scores_C);
assert(max(abs(source_culture_r_C))<1e-12); % Culture-orthogonal coincident sources
assert(rank(scores_M)==max_sources && rank(scores_C)==max_sources);
gram_M=(scores_M'*scores_M)/(nSubjects-1);
gram_C=(scores_C'*scores_C)/(nSubjects-1);
fprintf('Per-source culture correlations (M, C):\n');disp([source_culture_r_M;source_culture_r_C]);
X_M=[ones(nSubjects,1),vec_targ_real_M];
X_C=[ones(nSubjects,1),vec_targ_real_C];
tasks_ibin=repmat((1:nBins)',1,nSims_per_bin);tasks_ibin=tasks_ibin(:);
tasks_isim=repmat(1:nSims_per_bin,nBins,1);tasks_isim=tasks_isim(:);
rng(20260915,'twister');
scale_centers=randn(nCenters_M,nROIs);
scale_centers_M=scale_centers;
scale_centers_C=scale_centers-Q*(Q'*scale_centers);
dm=scale_centers_M(artifG_M,:); dm=dm-mean(dm,1);
dc=scale_centers_C(artifG_C,:); dc=dc-mean(dc,1);
delta_M=1+0.5*dm./max(abs(dm),[],1);
delta_C=1+0.5*dc./max(abs(dc),[],1);
scaling_culture_r_M=corr(vec_targ_real_M,delta_M);
scaling_culture_r_C=corr(vec_targ_real_C,delta_C);
assert(max(abs(scaling_culture_r_C))<1e-12); % Culture-orthogonal coincident sources
fprintf('Scaling max abs culture r M/C: %.6f %.3g\n',max(abs(scaling_culture_r_M)),max(abs(scaling_culture_r_C)));
assert(min([delta_M(:);delta_C(:)])>=0.5-1e-12 && max([delta_M(:);delta_C(:)])<=1.5+1e-12);
use_residual_scaling=true;
if ~use_residual_scaling; delta_M=ones(size(delta_M));delta_C=ones(size(delta_C));end
max_workers=Inf;
if strcmp(getenv('V3_SMOKE'),'1'),max_workers=0;
elseif isempty(gcp('nocreate')),parpool(min(8,feature('numcores')));end
for icase=1:numel(source_counts)
 K=source_counts(icase);
 all_gt=nan(nTotal,nROIs); all_sM=all_gt;all_eM=all_gt;all_sC=all_gt;all_eC=all_gt;
 cGT=nan(nTotal,1);sGTM=cGT;eGTM=cGT;sSAM=cGT;eSAM=cGT;
 sGTC=cGT;eGTC=cGT;sSAC=cGT;eSAC=cGT;
 energy=nan(nTotal,3); scales=nan(nTotal,2); spatial_rank=nan(nTotal,1);
 source_maps=nan(nTotal,nROIs,max_sources); actual_source_SA=nan(nTotal,max_sources);
 tic;
 parfor (itask=1:nTotal,max_workers)
  ibin=tasks_ibin(itask);isim=tasks_isim(itask);
  rng(1000000+ibin*1000+isim,'Threefry');
  if is_null_sims(ibin,isim);target_rho=-0.1+0.2*rand;else;target_rho=-1+2*rand;end
  spatial_noise=randn(nROIs,1);
  spatial_noise=zscore(spatial_noise-ideal_gradient*(ideal_gradient\spatial_noise));
  gt=max(mu_pos+amp_pos*(target_rho*zscore(ideal_gradient)+sqrt(1-target_rho^2)*spatial_noise),0.02);
  artif_noise=randn(nROIs,1);
  artif_noise=zscore(artif_noise-ideal_gradient*(ideal_gradient\artif_noise));
  shape=artif_rhos(ibin)*zscore(ideal_gradient)+sqrt(1-artif_rhos(ibin)^2)*artif_noise;
  B=zeros(nROIs,max_sources);
  B(:,1)=base_artif+range_artif*(shape-min(shape))/(max(shape)-min(shape));
  noise=randn(nSubjects,nROIs)*noise_level; % Measurement noise
  lo=max(0,artif_rhos(ibin)-0.1); hi=min(1,artif_rhos(ibin)+0.1);
  mean_r=lo+(hi-lo)*rand;
  rho_sources=mean_r*ones(1,max_sources);
  for pair=2:2:max_sources-1
   d=(2*rand-1)*min(1-mean_r,1+mean_r);
   rho_sources(pair:pair+1)=mean_r+[d,-d];
  end
  g=zscore(ideal_gradient);
  for j=1:max_sources
   u=randn(nROIs,1);u=u-mean(u);u=zscore(u-g*((g'*u)/(g'*g)));
   shape=rho_sources(j)*g+sqrt(max(0,1-rho_sources(j)^2))*u;
   B(:,j)=base_artif+range_artif*(shape-min(shape))/(max(shape)-min(shape));
  end
  assert(abs(mean(corr(ideal_gradient,B(:,1:K)))-mean_r)<1e-12);
  target_energy=mean(B(:,1).^2);
  BK=B(:,1:K);
  preM=mean(sum((BK*gram_M(1:K,1:K)).*BK,2));
  preC=mean(sum((BK*gram_C(1:K,1:K)).*BK,2));
  scaleM=sqrt(target_energy/preM);scaleC=sqrt(target_energy/preC);
  artifactM=scaleM*scores_M(:,1:K)*BK';
  artifactC=scaleC*scores_C(:,1:K)*BK';
  actualM=mean(sum(artifactM.^2,1)/(nSubjects-1));
  actualC=mean(sum(artifactC.^2,1)/(nSubjects-1));
  assert(max(abs([actualM actualC]/target_energy-1))<1e-10);
  YM=vec_targ_real_M*gt'+artifactM+delta_M.*noise;
  YC=vec_targ_real_C*gt'+artifactC+delta_C.*noise;
  bsM=zeros(nROIs,1);beM=bsM;bsC=bsM;beC=bsM;
  for r=1:nROIs
   y=YM(:,r);b=X_M\zscore(y);bsM(r)=b(2);
   w=(clean_pool_real_M'*y)/(nSubjects-1);
   [~,PC,~,~,exp_r]=pca(clean_pool_real_M.*w','NumComponents',max_pcs_M);
   XC=[ones(nSubjects,1),PC(:,1:min(max_pcs_M,length(exp_r)))];
   res=y-XC*(XC\y);b=X_M\zscore(res);beM(r)=b(2);
   y=YC(:,r);b=X_C\zscore(y);bsC(r)=b(2);
   w=(clean_pool_real_C'*y)/(nSubjects-1);
   [~,PC,~,~,exp_r]=pca(clean_pool_real_C.*w','NumComponents',max_pcs_C);
   XC=[ones(nSubjects,1),PC(:,1:min(max_pcs_C,length(exp_r)))];
   res=y-XC*(XC\y);b=X_C\zscore(res);beC(r)=b(2);
  end
  cGT(itask)=corr(gt,ideal_gradient);
  sGTM(itask)=corr(bsM,gt);eGTM(itask)=corr(beM,gt);
  sGTC(itask)=corr(bsC,gt);eGTC(itask)=corr(beC,gt);
  sSAM(itask)=corr(bsM,ideal_gradient);eSAM(itask)=corr(beM,ideal_gradient);
  sSAC(itask)=corr(bsC,ideal_gradient);eSAC(itask)=corr(beC,ideal_gradient);
  all_gt(itask,:)=gt';all_sM(itask,:)=bsM';all_eM(itask,:)=beM';
  all_sC(itask,:)=bsC';all_eC(itask,:)=beC';
  source_maps(itask,:,:)=reshape(B,[1 nROIs max_sources]);
  actual_source_SA(itask,:)=corr(ideal_gradient,B);
  energy(itask,:)=[target_energy actualM actualC];
  scales(itask,:)=[scaleM scaleC];spatial_rank(itask)=rank(BK);
 end
 corr_GT_SA=reshape(cGT,nBins,nSims_per_bin);
 corr_est_start_GT_M=reshape(sGTM,nBins,nSims_per_bin);
 corr_est_end_GT_M=reshape(eGTM,nBins,nSims_per_bin);
 corr_est_start_GT_C=reshape(sGTC,nBins,nSims_per_bin);
 corr_est_end_GT_C=reshape(eGTC,nBins,nSims_per_bin);
 corr_est_start_SA_M=reshape(sSAM,nBins,nSims_per_bin);
 corr_est_end_SA_M=reshape(eSAM,nBins,nSims_per_bin);
 corr_est_start_SA_C=reshape(sSAC,nBins,nSims_per_bin);
 corr_est_end_SA_C=reshape(eSAC,nBins,nSims_per_bin);
 maps_GT=all_gt;maps_sM=all_sM;maps_eM=all_eM;maps_sC=all_sC;maps_eC=all_eC;
 save(fullfile(out_dir,sprintf('multi_artifact_K%d.mat',K)),...
  'use_residual_scaling','delta_M','delta_C','scaling_culture_r_M','scaling_culture_r_C',...
  'K','source_counts','source_culture_r_M','source_culture_r_C','scores_M','scores_C',...
  'center_sources_raw','center_sources_M','center_sources_C','scale_centers_M','scale_centers_C',...
  'vec_targ_real_M','vec_targ_real_C','artifG_M','artifG_C','targG_M','targG_C','gram_M','gram_C','source_maps','actual_source_SA','energy','scales','spatial_rank',...
  'artif_rhos','nSims_per_bin','nNull_per_bin','is_null_sims','corr_GT_SA',...
  'nProvinces_M','nCenters_M','max_pcs_M','nProvinces_C','nCenters_C','max_pcs_C',...
  'orth_thres','clean_pool_real_M','clean_pool_real_C',...
  'corr_est_start_GT_M','corr_est_end_GT_M','corr_est_start_GT_C','corr_est_end_GT_C',...
  'corr_est_start_SA_M','corr_est_end_SA_M','corr_est_start_SA_C','corr_est_end_SA_C',...
  'maps_GT','maps_sM','maps_eM','maps_sC','maps_eC','-v7');
 fprintf('K=%d COMPLETED %.2f min | M %.4f -> %.4f | C %.4f -> %.4f | spatial rank [%d,%d]\n',...
  K,toc/60,mean(sGTM),mean(eGTM),mean(sGTC),mean(eGTC),min(spatial_rank),max(spatial_rank));
end
disp('MULTI_ARTIFACT_COMPLETED');





