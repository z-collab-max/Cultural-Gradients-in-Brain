function result=fit_china_mediation(T,cortex,scores,n_boot,seed,additional)
if nargin<4,n_boot=0;end
if nargin<5,seed=20261103;end
if nargin<6,additional={};end
rng(seed);R=size(cortex,2);tasks=numel(scores);
result.a=cell(2,tasks);result.b=cell(2,tasks);result.indirect=cell(2,tasks);
result.indirect_se=cell(2,tasks);result.indirect_ci=cell(2,tasks);result.total=cell(2,tasks);result.direct=cell(2,tasks);
result.behavior_effect=nan(R,tasks);result.task_indirect=nan(R,tasks,2);result.task_indirect_z=nan(R,tasks,2);
for t=1:tasks
 Y=scores{t};brain_beta=zeros(R,size(Y,2));
 for pc=1:2
  x=sprintf('culture_pc%d',pc);a=zeros(R,size(Y,2));b=a;c=a;direct=a;se=nan(size(a));ci=nan(R,size(Y,2),2);
  for r=1:R
   for q=1:size(Y,2)
    D=T;D.mediator=zscore(cortex(:,r));D.response=Y(:,q);
    required=[china_covariates,additional,{x,'mediator','response'}];D=D(all(isfinite(D{:,required}),2),:);
    ma=fitlme(D,china_formula('mediator',{x},additional),'FitMethod','REML');
    mb=fitlme(D,china_formula('response',{x,'mediator'},additional),'FitMethod','REML');
    mc=fitlme(D,china_formula('response',{x},additional),'FitMethod','REML');
    a(r,q)=estimate(ma,x);b(r,q)=estimate(mb,'mediator');c(r,q)=estimate(mc,x);direct(r,q)=estimate(mb,x);
    if pc==1
     assoc=fitlme(D,china_formula('response',{'mediator'},additional),'FitMethod','REML');
     brain_beta(r,q)=estimate(assoc,'mediator');
    end
    if n_boot>1
     draws=zeros(n_boot,1);groups=categories(removecats(D.province_birth));
     for k=1:n_boot
      chosen=randi(numel(groups),numel(groups),1);blocks=cell(numel(groups),1);
      for h=1:numel(groups)
       Q=D(D.province_birth==groups{chosen(h)},:);
       Q.province_city=categorical("B"+h+"_"+string(Q.province_city));Q.province_birth=categorical(repmat("B"+h,height(Q),1));blocks{h}=Q;
      end
      boot=vertcat(blocks{:});aa=fitlme(boot,china_formula('mediator',{x},additional),'FitMethod','REML');
      bb=fitlme(boot,china_formula('response',{x,'mediator'},additional),'FitMethod','REML');
      draws(k)=estimate(aa,x)*estimate(bb,'mediator');
     end
     se(r,q)=std(draws);ci(r,q,:)=quantile(draws,[.025 .975]);
    end
   end
  end
  result.a{pc,t}=a;result.b{pc,t}=b;result.total{pc,t}=c;result.direct{pc,t}=direct;result.indirect{pc,t}=a.*b;
  result.indirect_se{pc,t}=se;result.indirect_ci{pc,t}=ci;
  result.task_indirect(:,t,pc)=mean(a.*b,2);
  if n_boot>1,result.task_indirect_z(:,t,pc)=mean((a.*b)./se,2);end
 end
 result.behavior_effect(:,t)=mean(brain_beta,2);
 fprintf('Mediation task %d/%d complete\n',t,tasks);
end
result.bootstrap_repetitions=n_boot;
result.formulas={china_formula('mediator',{'culture_pc1'},additional),china_formula('response',{'culture_pc1','mediator'},additional),china_formula('response',{'culture_pc1'},additional)};
end
function value=estimate(fit,term)
j=find(strcmp(fit.Coefficients.Name,term));assert(isscalar(j));value=fit.Coefficients.Estimate(j);
end
