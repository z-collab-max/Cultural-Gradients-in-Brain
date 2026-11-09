function generate_china_resampling(config_file,output_dir)
% Refit the shared mixed model in each province or participant subsample.
cfg=jsondecode(fileread(config_file));if ~isfolder(output_dir),mkdir(output_dir);end
if ~isfield(cfg,'seed'),cfg.seed=20261029;end
if ~isfield(cfg,'replicates'),cfg.replicates=100;end
if ~isfield(cfg,'ratios'),cfg.ratios=0.1:0.1:0.9;end
rng(cfg.seed,'twister');S=load(cfg.brain_file,'brain_data');b=S.brain_data;
S=load(cfg.raw_file,'raw_res','valididx');assert(isfield(S,'valididx'),'Residual input requires valididx');
final_valid_idx=S.valididx(:);T=study_china_table(b,find(final_valid_idx));
Y=S.raw_res{1,29};assert(size(Y,1)==height(T));
C=[b.pc7_1(final_valid_idx) b.pc7_2(final_valid_idx)];
f=study_china_fit(T,Y,C);map_target_1=f.beta(1,:)';map_target_2=f.beta(2,:)';
prov_ids=findgroups(T.province_birth);unique_provs=unique(prov_ids);n_provs=numel(unique_provs);
ratios=cfg.ratios;n_perms=cfg.replicates;nr=numel(ratios);R=size(Y,2);
Beta_Prov_PC1=cell(nr,1);Beta_Prov_PC2=cell(nr,1);Beta_Subj_PC1=cell(nr,1);Beta_Subj_PC2=cell(nr,1);
Res_Prov_Corr1=nan(nr,n_perms);Res_Prov_Corr2=Res_Prov_Corr1;
Res_Subj_Corr1=Res_Prov_Corr1;Res_Subj_Corr2=Res_Prov_Corr1;
Draw_Prov=cell(nr,n_perms);Draw_Subj=Draw_Prov;fit_failures=struct('scheme',{},'ratio',{},'draw',{},'message',{});
for a=1:nr
 P1=nan(R,n_perms);P2=P1;S1=P1;S2=P1;
 for k=1:n_perms
  selected=randsample(unique_provs,min(n_provs,max(3,floor(n_provs*ratios(a)))));
  ip=find(ismember(prov_ids,selected));is=[];
  for g=unique_provs'
   rows=find(prov_ids==g);chosen=randsample(rows,max(1,floor(numel(rows)*ratios(a))));is=[is;chosen];
  end
  Draw_Prov{a,k}=ip;Draw_Subj{a,k}=is;
  for scheme=1:2
   if scheme==1,rows=ip;else,rows=is;end
   try
    f=study_china_fit(T(rows,:),Y(rows,:),C(rows,:));
    if scheme==1,P1(:,k)=f.beta(1,:)';P2(:,k)=f.beta(2,:)';else,S1(:,k)=f.beta(1,:)';S2(:,k)=f.beta(2,:)';end
   catch err
    fit_failures(end+1)=struct('scheme',scheme,'ratio',ratios(a),'draw',k,'message',err.message);
   end
  end
 end
 Beta_Prov_PC1{a}=P1;Beta_Prov_PC2{a}=P2;Beta_Subj_PC1{a}=S1;Beta_Subj_PC2{a}=S2;
 Res_Prov_Corr1(a,:)=corr(map_target_1,P1);Res_Prov_Corr2(a,:)=corr(map_target_2,P2);
 Res_Subj_Corr1(a,:)=corr(map_target_1,S1);Res_Subj_Corr2(a,:)=corr(map_target_2,S2);
end
save(fullfile(output_dir,'Robustness_Resampling_Betas_and_Corr.mat'),...
 'Res_Prov_Corr1','Res_Prov_Corr2','Res_Subj_Corr1','Res_Subj_Corr2',...
 'Beta_Prov_PC1','Beta_Prov_PC2','Beta_Subj_PC1','Beta_Subj_PC2','ratios','n_perms',...
 'map_target_1','map_target_2','Draw_Prov','Draw_Subj','fit_failures','-v7.3');
save(fullfile(output_dir,'generation_inputs.mat'),'cfg','final_valid_idx');
end
