function result=fit_china_mixed(T,cortex,additional)
if nargin<3,additional={};end
result=struct('beta',nan(size(cortex,2),2),'se',nan(size(cortex,2),2),'p',nan(size(cortex,2),2),'formulas',{{}});
for pc=1:2
 assert(ismember(size(cortex,3),[1,2]));Y=cortex(:,:,min(pc,size(cortex,3)));
 x=sprintf('culture_pc%d',pc);formula=china_formula('response',{x},additional);result.formulas{pc}=formula;
 for r=1:size(cortex,2)
  D=T;D.response=zscore(Y(:,r));fit=fitlme(D,formula,'FitMethod','REML');
  j=find(strcmp(fit.Coefficients.Name,x));result.beta(r,pc)=fit.Coefficients.Estimate(j);
  result.se(r,pc)=fit.Coefficients.SE(j);result.p(r,pc)=fit.Coefficients.pValue(j);
 end
end
end
