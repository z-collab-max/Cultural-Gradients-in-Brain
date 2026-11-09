function render_brain_icc(source_dir,out)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','figures'));
a=load(fullfile(source_dir,'ExtFig1_ICC_SourceData.mat'));b=load(fullfile(source_dir,'ICC_calculation.mat'),'calculation');
if ~isfolder(out),mkdir(out);end
for modality=["thickness","area","volume"]
    roi_names=regexprep(string(b.calculation.(modality).brain_var(:)),'_(area|volume)$','_thickness');
    for metric=["icc1","icc3"]
        values=max(0,a.ExtFig1_ICC_Data.(modality).(metric)(:));
        if metric=="icc1",color_limits=[-.01 .01];else,color_limits=[-1 1];end
        stem=char(modality+"_"+metric);file=fullfile(out,[stem '.mat']);save(file,'roi_names','values','color_limits');
        render_brain_map(file,fullfile(out,stem));
    end
end
end
