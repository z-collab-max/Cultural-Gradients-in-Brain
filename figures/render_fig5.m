function render_fig5(source_dir,out,draw_surface)
if nargin<3,draw_surface=true;end
addpath(fileparts(mfilename('fullpath')));
if ~isfolder(out),mkdir(out);end
models={'PC1','PC2','USA_TL'};letters='dchglk';
targets={'Eigenmode 2','Culture PC1 effect (China)','S–A axis','Culture PC2 effect (China)','Eigenmode 2','TL effect (China)'};
panels={};
for k=1:3
    model=models{k};brain=load(fullfile(source_dir,['Fig5_' model '_brain.mat']));
    for j=1:2
        i=2*k-2+j;p=load(fullfile(source_dir,sprintf('Fig5_%s_comparison%d.mat',model,j)));
        p.p_permutation=permutation_p_absolute(p.permutation_null,p.r);p.p_spin=permutation_p_absolute(p.spin_null,p.r);
        assert(isequal(string(p.roi(:)),string(brain.roi(:)))&&isequal(p.y(:),brain.beta(:)));
        panels{end+1}=struct('folder','Figure5','name',['Fig5' letters(i) '_scatter'],'kind','scatter',...
            'x',p.x(:),'y',p.y(:),'lh',startsWith(string(p.roi(:)),'lh_'),'xlabel',targets{i},'ylabel',[model ' effect'],...
            'stats',sprintf('r = %.3f, p_{spin} = %.4f',p.r,p.p_spin));
        panels{end+1}=struct('folder','Figure5','name',['Fig5' letters(i) '_permutation'],'kind','histogram',...
            'null',p.permutation_null(:),'observed',p.r,'p',p.p_permutation,'xlabel','Spatial correlation');
    end
    if draw_surface
        roi_names=brain.roi;values=brain.t;limit=10;if k==3,limit=14;end;color_limits=[-limit limit];
        file=fullfile(source_dir,['Fig5_' model '_surface.mat']);save(file,'roi_names','values','color_limits','-v7');
        render_brain_map(file,fullfile(out,['Fig5_' model '_brain']));
    end
end
file=fullfile(source_dir,'Fig5_panels.mat');save(file,'panels','-v7');render_panels(file,out);
end
