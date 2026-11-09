function render_rice_models(source_dir,output_dir)
figures=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'figures');
addpath(figures);
if ~isfolder(output_dir),mkdir(output_dir);end
for phenotype=["thickness","area"]
    for indicator=["PercentPaddy","RiceSuitability"]
        key=phenotype+"_"+indicator;
        d=load(fullfile(source_dir,"EDFig8_"+key+".mat"));
        for test=["spin","province"]
            panel=struct('folder','ExtendedDataFigure8','name',char(key+"_"+test),'kind','violin',...
                'beta',d.observed,'null',d.(test+"_null"),'p',d.("p_"+test),'labels',{d.gradient_names});
            panels={panel};
            input=fullfile(output_dir,key+"_"+test+"_panel.mat");
            save(input,'panels');render_panels(input,output_dir);
        end
        map=struct('roi_names',{d.roi_names},'values',d.tVal,'significant',d.fdr_significant,...
            'color_limits',[-1 1]*(4+(max(abs(d.tVal))>4)));
        input=fullfile(output_dir,key+"_surface.mat");save(input,'-struct','map');
        render_brain_map(input,char(fullfile(output_dir,'ExtendedDataFigure8',key+"_surface")));
    end
end
end


