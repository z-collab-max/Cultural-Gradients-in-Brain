function fig2_recovery_mixed(input_file,output_file)
% City-stratified residual bootstrap with province/city random intercepts.
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'study_inputs'),fullfile(here,'..','china'));
S=load(input_file);Y=double(S.Y);N=size(Y,1);R=size(Y,2);B=double(S.draws);
assert(size(S.covariates,2)==28 && isequal(size(S.scores),[N 2]) && B>=3);
assert(all(isfinite([Y S.covariates S.scores S.groups]),'all'));
names=china_covariates;b=struct;for j=1:28,b.(names{j})=S.covariates(:,j);end
b.province_birth=S.groups(:,1);b.province_city=S.groups(:,2);T=study_china_table(b);
gp=findgroups(T.province_birth);gc=findgroups(T.province_city);P=max(gp);C=max(gc);
assert(P>=2 && C>P,'Two random intercepts require multiple provinces and cities within provinces');
counts=accumarray(gc,1);assert(all(counts>=2),'Each city needs at least two observations for residual resampling');
Zp=sparse((1:N)',gp,1,N,P);Zc=sparse((1:N)',gc,1,N,C);
formula=china_formula('response',{'culture'});calibration_errors=zeros(B,R,2);validation_errors=calibration_errors;
observed_beta=zeros(2,R);observed_se=observed_beta;residual_variance=observed_beta;
province_variance=observed_beta;city_variance=observed_beta;
calibration_seed=[20260915 20260916];validation_seed=[2026091601 2026091602];
for pc=1:2
 D=T;D.culture=zscore(S.scores(:,pc));
 assert(rank([ones(N,1) D{:,names} D.culture])==30,'Rank-deficient fixed design');
 nuisance=zeros(N,R);pool=zeros(N,R);up=zeros(P,R);uc=zeros(C,R);
 for r=1:R
  D.response=Y(:,r);m=fitlme(D,formula,'FitMethod','REML');
  X=designMatrix(m,'Fixed');coefs=m.Coefficients.Estimate;ix=strcmp(m.Coefficients.Name,'culture');
  observed_beta(pc,r)=coefs(ix);observed_se(pc,r)=m.Coefficients.SE(ix);
  [Psi,sigma2]=covarianceParameters(m);groups=m.Formula.GroupingVariableNames;
  vp=[];vc=[];
  for k=1:numel(groups)
   label=string(groups{k});assert(isscalar(label));
   if label=="province_birth",vp=Psi{k};elseif label=="province_city",vc=Psi{k};end
  end
  assert(isscalar(vp)&&isscalar(vc)&&sigma2>0);
  residual_variance(pc,r)=sigma2;province_variance(pc,r)=vp;city_variance(pc,r)=vc;
  fixed=X*coefs;V=sigma2*eye(N)+vp*full(Zp*Zp')+vc*full(Zc*Zc');weighted=V\(Y(:,r)-fixed);
  up(:,r)=vp*(Zp'*weighted);uc(:,r)=vc*(Zc'*weighted);
  nuisance(:,r)=fixed-D.culture*coefs(ix);
  e=Y(:,r)-fixed-Zp*up(:,r)-Zc*uc(:,r);e=e-Zc*(Zc'*e./counts);
  assert(mean(e.^2)>0);pool(:,r)=e*sqrt(sigma2/mean(e.^2));
 end
 Fp=random_factor(up,province_variance(pc,:));Fc=random_factor(uc,city_variance(pc,:));
 for batch=1:2
  if batch==1,rng(calibration_seed(pc),'twister');else,rng(validation_seed(pc),'twister');end
  errors=zeros(B,R);
  for draw=1:B
   idx=zeros(N,1);
   for g=1:C,rows=find(gc==g);idx(rows)=rows(randi(numel(rows),numel(rows),1));end
   pdraw=randn(P,R)*Fp;cdraw=randn(C,R)*Fc;
   simulated=nuisance+pdraw(gp,:)+cdraw(gc,:)+pool(idx,:);
   for r=1:R
    D.response=simulated(:,r);m=fitlme(D,formula,'FitMethod','REML');
    ix=strcmp(m.Coefficients.Name,'culture');errors(draw,r)=m.Coefficients.Estimate(ix);
   end
  end
  if batch==1,calibration_errors(:,:,pc)=errors;else,validation_errors(:,:,pc)=errors;end
 end
 fprintf('PC%d REML recovery: %d ROIs, %d draws per batch\n',pc,R,B);
end
assert(all(isfinite([calibration_errors(:);validation_errors(:)])));
save(output_file,'calibration_errors','validation_errors','observed_beta','observed_se',...
 'residual_variance','province_variance','city_variance','calibration_seed','validation_seed','formula','-v7');
end

function F=random_factor(u,variance)
% Shrink the cross-ROI BLUP correlation, preserving fitted marginal variances.
u=u-mean(u,1);scale=sqrt(sum(u.^2,1));v=zeros(size(u));ok=scale>1e-12;v(:,ok)=u(:,ok)./scale(ok);
C=v'*v;C(1:size(C,1)+1:end)=1;C=.5*C+.5*eye(size(C));
[Q,L]=eig((C+C')/2);F=diag(sqrt(max(diag(L),0)))*Q'*diag(sqrt(max(variance,0)));
assert(max(abs(sum(F.^2,1)-variance))<1e-7*max(1,max(variance)));
end
