function run_combat_comparison(input_dir,out,outer_draws,inner_draws,seed)
% Fit both protection designs with identical eligible rows and shuffle indices.
if nargin<3,outer_draws=10000;end
if nargin<4,inner_draws=10000;end
if nargin<5,seed=42;end
for protect=[false,true]
    if protect,name='protect_PC12_recruit';else,name='protect_PC12';end
    run_china_chain(input_dir,fullfile(out,name),outer_draws,inner_draws,seed,protect);
end
a=load(fullfile(out,'protect_PC12','main_maps.mat'));
b=load(fullfile(out,'protect_PC12_recruit','main_maps.mat'));
assert(isequal(a.roi_names,b.roi_names));
roi_names=a.roi_names;beta_pc12=a.beta;beta_recruit=b.beta;
null_pc12=a.null_beta;null_recruit=b.null_beta;
r=arrayfun(@(pc)corr(a.beta(pc,:)',b.beta(pc,:)'),1:2);
save(fullfile(out,'comparison.mat'),'roi_names','beta_pc12','beta_recruit','null_pc12','null_recruit','r','-v7');
panels=cell(1,2);
for pc=1:2
    panels{pc}=struct('folder','ExtendedDataFigure5','name',sprintf('Protection_PC%d',pc),...
        'kind','combat_comparison','x',a.beta(pc,:),'y',b.beta(pc,:),...
        'xlabel','Protect PC1 and PC2','ylabel','Also protect recruitment scores');
end
file=fullfile(out,'panels.mat');save(file,'panels');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','figures'));render_panels(file,out);
end
