function run_model_branches(input_file,output_dir)
% Input: branches struct array with name, T, cortex, predictor_names.
S=load(input_file);assert(isfield(S,'branches'));
if ~isfolder(output_dir),mkdir(output_dir);end
names={S.branches.name};assert(numel(unique(names))==numel(names));
for k=1:numel(S.branches)
 b=S.branches(k);assert(~isempty(regexp(b.name,'^[A-Za-z][A-Za-z0-9_]*$','once')));
 T=prepare_china_table(b.T);predictors=cellstr(string(b.predictor_names));M=numel(predictors);R=size(b.cortex,2);
 beta=zeros(M,R);p=beta;se=beta;extra={};if isfield(b,'additional_covariates'),extra=cellstr(string(b.additional_covariates));end
 for j=1:M
  T.(predictors{j})=zscore(double(T.(predictors{j})));
  for r=1:R
   T.response=zscore(double(b.cortex(:,r)));m=fitlme(T,china_formula('response',predictors(j),extra),'FitMethod','REML');
   ix=find(strcmp(m.Coefficients.Name,predictors{j}));assert(isscalar(ix));beta(j,r)=m.Coefficients.Estimate(ix);se(j,r)=m.Coefficients.SE(ix);p(j,r)=m.Coefficients.pValue(ix);
  end
 end
 model_labels=cellstr(string(b.name)+"_"+string(predictors));save(fullfile(output_dir,[b.name '.mat']),'beta','p','se','model_labels','-v7');
end
end
