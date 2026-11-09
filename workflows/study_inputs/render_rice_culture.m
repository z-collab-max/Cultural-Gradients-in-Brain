function render_rice_culture(input,out)
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'..','..','figures'));
d=jsondecode(fileread(input));rows=d.rows;
if ~iscell(rows{1}), rows=mat2cell(rows,ones(size(rows,1),1),size(rows,2));end
labels=unique(string(cellfun(@(r) r{2},rows,'UniformOutput',false)),'stable');
indicators=unique(string(cellfun(@(r) r{1},rows,'UniformOutput',false)),'stable');
r=nan(numel(labels),numel(indicators));p=r;
for j=1:numel(rows)
    row=rows{j};i=labels==row{2};k=indicators==row{1};r(i,k)=row{3};p(i,k)=row{4};
end
d=struct('folder','ExtendedDataFigure8','name','Rice_culture','kind','heatmap','matrix',r,'p',p,'labels',labels,'xlabels',indicators,'clim',[-1 1]);
panels={d};if ~isfolder(out),mkdir(out);end
file=fullfile(out,'rice_culture_panel.mat');save(file,'panels');render_panels(file,out);
end

