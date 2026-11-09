function generate_china_trajectory_inputs(brain_file,province_file,shuffle_file,out)
% Center-shuffle artifact components and residuals at each removal count.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
if ~isfolder(out),mkdir(out);end
S=load(brain_file,'brain_data');b=S.brain_data;
tp=readtable(province_file,Sheet='Sheet2');cul=table2struct(tp,'ToScalar',true);
covs=china_covariates;
allvars=fieldnames(b);brain_var=allvars(contains(allvars,'thickness') & ~contains(allvars,{'500_aparc_','unknown','WhiteSurfArea'}));
Y=cols(b,brain_var);X=cols(b,covs);birth=cols(b,{'pc7_1','pc7_2'});
extra=cols(b,{'centerID_rerange','province_birth','province_city'});
valididx=all(isfinite([Y X birth extra]),2);Y=Y(valididx,:);X=X(valididx,:);birth=birth(valididx,:);MRIID=b.MRIID(valididx);
centers=b.centerID_rerange(valididx);[~,first,group]=unique(centers);counts=accumarray(group,1);n=numel(group);
sh=load(shuffle_file,'center_ids','culture_labels');[found,center_group]=ismember(centers,sh.center_ids);assert(all(found));
assert(size(sh.culture_labels,2)>=10000);raw_res=cell(2,29);design=[ones(n,1),zscore(X)];Y=zscore(Y);
for pc=1:2
 name=sprintf('pc7_%d',pc);raw_scores=cul.(name)(sh.culture_labels(center_group,1:10000));scores=zscore(raw_scores);
 rho=corr(scores,birth(:,pc));varying=(max(raw_scores,[],1)>min(raw_scores,[],1))';
 retained=varying & abs(rho)<.01;scores=scores(:,retained);
 nuisance=design(:,2:end);coef=zeros(size(Y,2),size(scores,2));
 for j=1:size(scores,2),A=[ones(n,1),scores(:,j),nuisance];v=A\Y;coef(:,j)=v(2,:)';end
 for k=1:29,raw_res{pc,k}=zeros(size(Y));end
 raw_res{pc,1}=Y;explained=zeros(size(Y,2),30);
 for r=1:size(Y,2)
  features=scores(first,:).*coef(r,:);features=features-counts'*features/n;
  [U,S,~]=svd(sqrt(counts).*features,'econ');sv=diag(S);rank_pc=sum(sv>max(size(features))*eps(max(sv)));
  basis=U(group,1:rank_pc)./sqrt(counts(group));contributions=basis.*(basis'*Y(:,r))';removed=cumsum(contributions,2);
  for k=2:29,raw_res{pc,k}(:,r)=zscore(Y(:,r)-removed(:,min(k,rank_pc)));end
  explained(r,1:numel(sv))=100*(sv.^2)'/sum(sv.^2);
 end
 save(fullfile(out,sprintf('artifact_PC%d.mat',pc)),'rho','retained','coef','explained','covs');
 fprintf('PC%d: %d individuals, %d retained shuffles, 29 residual levels\n',pc,n,sum(retained));
end
removed_pcs=[0,2:29];
save(fullfile(out,'raw_res_birth_all_expcs.mat'),'raw_res','valididx','MRIID','brain_var','removed_pcs','covs','-v7.3');
end

function M=cols(b,names)
M=zeros(numel(b.age),numel(names));for j=1:numel(names),M(:,j)=b.(names{j})(:);end
end
