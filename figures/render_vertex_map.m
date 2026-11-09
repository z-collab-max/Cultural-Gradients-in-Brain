function render_vertex_map(input_file,output_stem)
% MAT input: values_lh and values_rh in fsaverage5 vertex order, color_limits.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'simulation'));
d=load(input_file);limits=double(d.color_limits);assert(numel(limits)==2&&diff(limits)>0);
rgb=[33 102 172;146 197 222;247 247 247;244 165 130;178 24 43]/255;
cm=interp1(linspace(0,1,5),rgb,linspace(0,1,256));
f=figure('Visible','off','Color','w','Position',[100 100 1200 650]);
for h=1:2
 if h==1,hemi='lh';angles=[-90 90];else,hemi='rh';angles=[90 -90];end
 s=load(fullfile(root,'assets','surfaces',['fsaverage5_' hemi '.mat']));v=double(d.(['values_' hemi])(:));
 assert(numel(v)==size(s.vertices,1)&&any(isfinite(v(logical(s.cortex)))));
 q=1+round(255*min(1,max(0,(v-limits(1))/diff(limits))));q(~isfinite(q))=1;
 c=cm(q,:);mask=~logical(s.cortex(:)) | ~isfinite(v);c(mask,:)=repmat([.9 .9 .9],sum(mask),1);
 for side=1:2
  ax=axes(f,'Position',[.03+(side-1)*.49 .55-(h-1)*.46 .44 .40]);
  patch(ax,'Vertices',double(s.vertices),'Faces',double(s.faces),'FaceVertexCData',c,'FaceColor','interp','EdgeColor','none');
  axis(ax,'tight');axis(ax,'equal');axis(ax,'off');view(ax,angles(side),0);camtarget(ax,mean(double(s.vertices),1));camproj(ax,'orthographic');camlight(ax,'headlight');lighting(ax,'gouraud');material(ax,'dull');colormap(ax,cm);clim(ax,limits);
 end
end
cb=colorbar(ax,'Location','southoutside');cb.Position=[.3 .035 .4 .025];publication_style(f);
parent=fileparts(output_stem);if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
savefig(f,[output_stem '.fig']);exportgraphics(f,[output_stem '.png'],'Resolution',180);close(f);
end

