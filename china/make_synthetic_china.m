function [T,cortex,behavior,task_names]=make_synthetic_china(seed,n_rois)
if nargin<1,seed=20261103;end
if nargin<2,n_rois=3;end
rng(seed);n_provinces=20;cities_per_province=3;per_city=12;
city=repelem((1:n_provinces*cities_per_province)',per_city);
province=ceil(city/cities_per_province);n=numel(city);
T=table(province,city,'VariableNames',{'province_birth','city_birth'});
names=china_covariates;
for j=1:numel(names)
 if j<=17,T.(names{j})=randn(n,1);
 elseif j<=26,v=randn(max(city),1);T.(names{j})=v(city);
 else,v=randn(n_provinces,1);T.(names{j})=v(province);end
end
T.age=18+42*rand(n,1);T.sex=double(rand(n,1)>.5);
T.education=10+8*rand(n,1);T.father_education=4+14*rand(n,1);
T.mother_education=4+14*rand(n,1);
T.intracranial_volume=1500+130*randn(n,1);T.total_surface_area=180000+16000*randn(n,1);
scores=zscore(randn(n_provinces,2));T.culture_pc1=scores(province,1);T.culture_pc2=scores(province,2);
T.relocator=double(rand(n,1)>.5);destination=province;
moving=find(T.relocator==1);
destination(moving)=mod(province(moving)-1+randi(n_provinces-1,numel(moving),1),n_provinces)+1;
T.province_destination=destination;
T.destination_pc1=scores(destination,1);T.destination_pc2=scores(destination,2);
X=zscore(T{:,names});cortex=zeros(n,n_rois);
for r=1:n_rois
 u=.18*randn(n_provinces,1);v=.12*randn(max(city),1);
 cortex(:,r)=.2*sin(r)*T.culture_pc1+.18*cos(r)*T.culture_pc2+X*(.025*randn(28,1))+...
 u(province)+v(city)+.1*T.relocator.*T.culture_pc1+.08*(T.destination_pc1-T.culture_pc1)+randn(n,1);
end
task_names={'CVLT','GoNogo','NBack','Passball','RO','SDMT','UG'};
behavior=cell(1,7);
for t=1:7
 latent=.15*T.culture_pc1+.07*T.culture_pc2+(.08*t)*cortex(:,1)+X*(.02*randn(28,1))+randn(n,1);
 behavior{t}=latent*(.6+.4*rand(1,3))+.6*randn(n,3);
end
end
