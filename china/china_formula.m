function formula=china_formula(response,predictors,additional)
if nargin<3,additional={};end
terms=[predictors,china_covariates,additional];
formula=sprintf('%s ~ 1 + %s + (1|province_birth) + (1|province_city)',...
 response,strjoin(terms,' + '));
end
