function render_fig4_spatial(source_dir,out,draw_surface)
if nargin<3,draw_surface=false;end
for label={'behavior','PC1','PC2'}
    branch=label{1};folder=fullfile(source_dir,branch);render_panels(fullfile(folder,'panels.mat'),out);
    if draw_surface
        for figure={'fig4','ed12'}
            tag=figure{1};s=load(fullfile(folder,[tag '_source.mat']));roi_names=s.roi_names;
            if strcmp(tag,'fig4'),values=s.surface_z;color_limits=[-4 4];
            else,values=s.surface_slope;color_limits=s.color_limits;end
            file=fullfile(folder,[tag '_surface.mat']);save(file,'roi_names','values','color_limits','-v7');
            render_brain_map(file,fullfile(out,[tag '_' branch '_surface']));
        end
    end
end
end
