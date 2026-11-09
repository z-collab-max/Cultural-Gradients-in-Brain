function result=fit_china_relocation(T,cortex,additional)
if nargin<3,additional={};end
result=struct('beta',nan(size(cortex,2),4,2),'group_beta',nan(size(cortex,2),2,2),'formulas',{{}});
for pc=1:2
 assert(ismember(size(cortex,3),[1,2]));Y=cortex(:,:,min(pc,size(cortex,3)));
 x=sprintf('culture_pc%d',pc);delta=sprintf('culture_delta%d',pc);interaction=[x ':relocator'];
 formula=china_formula('response',{x,'relocator',interaction,delta},additional);result.formulas{pc}=formula;
 for r=1:size(cortex,2)
  D=T;D.response=zscore(Y(:,r));fit=fitlme(D,formula,'FitMethod','REML');
  terms={x,'relocator',interaction,delta};
  for k=1:4
   j=find(strcmp(fit.Coefficients.Name,terms{k}));
   if isempty(j)&&k==3,j=find(strcmp(fit.Coefficients.Name,['relocator:' x]));end
   assert(isscalar(j));result.beta(r,k,pc)=fit.Coefficients.Estimate(j);
  end
  for g=0:1
   sub=D(D.relocator==g,:);m=fitlme(sub,china_formula('response',{x},additional),'FitMethod','REML');
   j=find(strcmp(m.Coefficients.Name,x));result.group_beta(r,g+1,pc)=m.Coefficients.Estimate(j);
  end
 end
end
result.effect_names={'birthplace','relocation_shift','deculturation','acculturation'};
result.difference_map=squeeze(result.group_beta(:,1,:)-result.group_beta(:,2,:));
end
