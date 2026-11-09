function render_culture_rsa(source_dir,out)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','figures'));
s=load(fullfile(source_dir,'RSA_source_results.mat'));q=load(fullfile(source_dir,'generation_inputs.mat'),'brain_var');roi_names=q.brain_var;
panels={};if ~isfolder(out),mkdir(out);end
for pc=1:2
    d=s.ED5_bc_Data.(sprintf('PC%d',pc));stem=sprintf('RSA_PC%d',pc);
    values=d.brainmap_z;color_limits=[-1 1]*(7-2*pc);file=fullfile(out,[stem '_surface.mat']);
    save(file,'roi_names','values','color_limits');render_brain_map(file,fullfile(out,[stem '_surface']));
    r=corr(d.scatter_x(:),d.scatter_y(:));
    panels{end+1}=struct('folder','ExtendedDataFigure6','name',[stem '_culture'],'kind','scatter',...
        'x',d.scatter_x,'y',d.scatter_y,'lh',startsWith(string(roi_names),'lh_'),'xlabel','RSA','ylabel','Cultural effect',...
        'stat_null',d.spin_null,'stat_observed',r,'stat_label','r','stats',sprintf('r = %.3f, p_{spin} = %.3f',r,d.spin_p));
    p=mean(abs(d.hist_null_dist)>=abs(d.hist_real_val));
    panels{end+1}=struct('folder','ExtendedDataFigure6','name',[stem '_reference'],'kind','histogram',...
        'null',d.hist_null_dist,'observed',d.hist_real_val,'p',p,'xlabel',d.hist_ref_name);
end
file=fullfile(out,'rsa_panels.mat');save(file,'panels');render_panels(file,out);
end
