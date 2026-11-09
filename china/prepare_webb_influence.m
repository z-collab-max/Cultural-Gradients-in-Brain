function result=prepare_webb_influence(T,cortex,pc)
% Province-aggregated HC2 scores from the fitted two-intercept covariance.
T=prepare_china_table(T);name=sprintf('culture_pc%d',pc);
province=findgroups(T.province_birth);G=max(province);R=size(cortex,2);N=height(T);
result.point=zeros(R,1);result.influence=zeros(G,R);result.max_leverage=zeros(R,1);
for r=1:R
 D=T;D.response=zscore(cortex(:,r));m=fitlme(D,china_formula('response',{name}),'FitMethod','REML');
 X=designMatrix(m,'Fixed');[Psi,sigma2]=covarianceParameters(m);
 grouping=m.Formula.GroupingVariableNames;V=sigma2*eye(N);
 for k=1:numel(grouping)
  group_name=string(grouping{k});assert(isscalar(group_name));
  ids=findgroups(D.(group_name));Z=sparse((1:N)',ids,1,N,max(ids));
  assert(isscalar(Psi{k}));V=V+Psi{k}*full(Z*Z');
 end
 L=chol(V,'lower');Xw=L\X;yw=L\D.response;[Q,RR]=qr(Xw,0);
 b=RR\(Q'*yw);assert(max(abs(b-m.Coefficients.Estimate))<1e-5);
 h=sum(Q.^2,2);assert(max(h)<1);e=(yw-Xw*b)./sqrt(max(1-h,1e-7));
 slope=RR\Q';j=find(strcmp(m.Coefficients.Name,name));assert(isscalar(j));
 result.point(r)=b(j);score=slope(j,:)'.*e;
 result.influence(:,r)=accumarray(province,score,[G,1]);result.max_leverage(r)=max(h);
end
result.formula=china_formula('response',{name});result.adjustment='whitened_HC2';
end
