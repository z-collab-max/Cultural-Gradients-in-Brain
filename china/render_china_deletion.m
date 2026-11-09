function render_china_deletion(source,out)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','figures'));
d=load(source);kind=string(d.kind(:));labels=string(d.labels(:));status=string(d.status(:));
assert(numel(kind)==size(d.correlation,1));panels={};
for pc=1:2
    ix=kind=="LOO" & status=="ok";
    panels{end+1}=struct('folder','ExtendedDataFigure10','name',sprintf('China_PC%d_LOO',pc),...
        'kind','bars','y',d.correlation(ix,pc),'labels',{cellstr(labels(ix))},...
        'color',repmat([.44 .51 .62],sum(ix),1),'ylabel','Spatial correlation');
    ix=kind=="progressive" & status=="ok";
    if any(ix)
        panels{end+1}=struct('folder','ExtendedDataFigure10','name',sprintf('China_PC%d_deletion',pc),...
            'kind','curves','x',str2double(labels(ix)),'y',d.correlation(ix,pc),...
            'labels',{{'Progressive deletion'}},'xlabel','Provinces omitted','ylabel','Spatial correlation');
    end
end
if ~isfolder(out),mkdir(out);end
file=fullfile(out,'panels.mat');save(file,'panels');render_panels(file,out);
end
