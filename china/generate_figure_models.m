function result=generate_figure_models(input_file,output_file)
% Input MAT contains T, cortex and predictor_names; rows must be aligned.
S=load(input_file);T=prepare_china_table(S.T);Y=double(S.cortex);
assert(height(T)==size(Y,1)&&all(isfinite(Y),'all'));
names=cellstr(string(S.predictor_names));R=size(Y,2);M=numel(names);
result.beta=zeros(M,R);result.p=zeros(M,R);result.low=zeros(M,R);result.high=zeros(M,R);
for j=1:M
 assert(ismember(names{j},T.Properties.VariableNames),'Predictor is missing');
 assert(isnumeric(T.(names{j}))&&all(isfinite(T.(names{j}))));
 T.(names{j})=zscore(T.(names{j}));
 for r=1:R
  T.response=zscore(Y(:,r));m=fitlme(T,china_formula('response',names(j)),'FitMethod','REML');
  k=find(strcmp(m.Coefficients.Name,names{j}));assert(isscalar(k));ci=coefCI(m);
  result.beta(j,r)=m.Coefficients.Estimate(k);result.p(j,r)=m.Coefficients.pValue(k);
  result.low(j,r)=ci(k,1);result.high(j,r)=ci(k,2);
 end
end
result.predictor_names=names;result.n=height(T);
parent=fileparts(output_file);if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
save(output_file,'-struct','result');
end
