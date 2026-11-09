function coefficients=null_effect_coefficients(T,cortex,score_pool,pc)
% score_pool: individuals-by-shuffles, already screened for artifact estimation.
T=prepare_china_table(T);name=sprintf('culture_pc%d',pc);
assert(size(score_pool,1)==height(T));coefficients=zeros(size(score_pool,2),size(cortex,2));
for b=1:size(score_pool,2)
 D=T;D.(name)=zscore(score_pool(:,b));
 for r=1:size(cortex,2)
  D.response=zscore(cortex(:,r));fit=fitlme(D,china_formula('response',{name}),'FitMethod','REML');
  j=strcmp(fit.Coefficients.Name,name);coefficients(b,r)=fit.Coefficients.Estimate(j);
 end
end
end
