function result=study_china_fit(T,Y,predictors)
% REML variance components are estimated separately for each response and draw.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
assert(size(predictors,1)==height(T) && size(Y,1)==height(T));
assert(all(isfinite([Y predictors]),'all'));
R=size(Y,2);B=size(predictors,2);
for key={'beta','t','p','df'},result.(key{1})=nan(B,R);end
formula=china_formula('response',{'study_predictor'});
for k=1:B
 D=T;D.study_predictor=zscore(predictors(:,k));
 assert(rank([ones(height(D),1),D{:,china_covariates},D.study_predictor])==30,'Fixed effects are rank deficient');
 for r=1:R
  D.response=zscore(Y(:,r));m=fitlme(D,formula,'FitMethod','REML');
  c=m.Coefficients(strcmp(m.Coefficients.Name,'study_predictor'),:);
  result.beta(k,r)=c.Estimate;result.t(k,r)=c.tStat;
  result.p(k,r)=c.pValue;result.df(k,r)=c.DF;
 end
end
result.formula=formula;
end
