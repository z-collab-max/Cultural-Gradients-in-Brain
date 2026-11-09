function render_china_resampling(source,out)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','figures'));
s=load(source);panels=cell(1,2);
for pc=1:2
    data=cat(3,s.(sprintf('Prov_PC%d',pc)),s.(sprintf('Subj_PC%d',pc)));
    panels{pc}=struct('folder','ExtendedDataFigure9','name',sprintf('China_PC%d',pc),'kind','ribbons',...
        'x',s.ratios,'series',permute(data,[3 1 2]),'labels',{{'Province','Subject'}});
end
if ~isfolder(out),mkdir(out);end
file=fullfile(out,'resampling_panels.mat');save(file,'panels');render_panels(file,out);
end
