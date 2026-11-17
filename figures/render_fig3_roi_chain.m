function render_fig3_roi_chain(source_dir,out,draw_surface)
if nargin<3,draw_surface=false;end
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'layouts'));
plot_fig3jk_ed11b(fullfile(source_dir,'roi_source.mat'),out);
render_panels(fullfile(source_dir,'network_source.mat'),out);
if draw_surface
    d=load(fullfile(source_dir,'roi_source.mat'));roi_names=d.names;
    for pc=1:2
        values=d.roi_effect(:,pc).*(d.posidx(:,pc)|d.negidx(:,pc));color_limits=[-.15 .15];
        name=sprintf('Fig3i_PC%d',pc);file=fullfile(source_dir,[name '_surface.mat']);save(file,'roi_names','values','color_limits','-v7');
        render_brain_map(file,fullfile(out,name));
    end
end
end
