function formula=china_formula(response,predictors,additional,T)
if nargin<3,additional={};end
if nargin<4,controls=china_covariates;else,controls=china_covariates(T);end
terms=[predictors,controls,additional];
formula=sprintf('%s ~ 1 + %s + (1|province_birth) + (1|province_city)',...
 response,strjoin(terms,' + '));
end
