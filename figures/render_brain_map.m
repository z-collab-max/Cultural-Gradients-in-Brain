function render_brain_map(input_file,output_stem)
% MAT input: roi_names, values, color_limits; optional significant ROI mask.
here=fileparts(mfilename('fullpath'));root=fileparts(here);addpath(fullfile(root,'simulation'));addpath(fullfile(here,'helpers'));
d=load(input_file);names=string(d.roi_names(:));values=double(d.values(:));
d.color_limits=double(d.color_limits);
assert(numel(names)==numel(values)&&numel(unique(names))==numel(names));
assert(all(isfinite(values))&&numel(d.color_limits)==2&&diff(d.color_limits)>0);
rgb=[33 102 172;146 197 222;247 247 247;244 165 130;178 24 43]/255;
cm=interp1(linspace(0,1,5),rgb,linspace(0,1,256));
f=figure('Visible','off','Color','w','Position',[100 100 1200 650]);found=false(size(names));
for h=1:2
 if h==1,hemi='lh';cap='L';angles=[-90 90];else,hemi='rh';cap='R';angles=[90 -90];end
 s=load(fullfile(root,'assets','surfaces',['atlas_' cap '.mat']));
 label_names=string(s.label_names(:));ids=double(s.label_ids(:));v=nan(numel(s.label),1);sig=[];
 for r=1:numel(names)
  if ~startsWith(names(r),[hemi '_']),continue;end
  short=regexprep(names(r),'^[lr]h_|_(thickness|area|volume)$','');j=find(label_names==short);
  assert(isscalar(j),'ROI does not match exactly one atlas label: %s',names(r));
  v(double(s.label)==ids(j))=values(r);found(r)=true;
  if isfield(d,'significant')&&d.significant(r),sig(end+1)=ids(j);end
 end
 q=1+round(255*min(1,max(0,(v-d.color_limits(1))/diff(d.color_limits))));q(~isfinite(q))=1;
 colors=cm(q,:);colors(~isfinite(v),:)=repmat([.9 .9 .9],sum(~isfinite(v)),1);
 for side=1:2
  ax=axes(f,'Position',[.03+(side-1)*.49 .55-(h-1)*.46 .44 .40]);
  patch(ax,'Vertices',s.vertices,'Faces',double(s.faces),'FaceVertexCData',colors,'FaceColor','interp','EdgeColor','none');hold(ax,'on');
  for id=sig
   contours=roi_midpoint_segments(s.vertices,double(s.faces),double(s.label)==id);
   for b=1:numel(contours),vtx=contours{b};plot3(ax,vtx(:,1),vtx(:,2),vtx(:,3),'Color',[.25 .3 .35],'LineWidth',.85);end
  end
  axis(ax,'tight');axis(ax,'equal');axis(ax,'off');view(ax,angles(side),0);camtarget(ax,mean(s.vertices,1));camproj(ax,'orthographic');camlight(ax,'headlight');lighting(ax,'gouraud');material(ax,'dull');
  colormap(ax,cm);clim(ax,d.color_limits);
 end
end
assert(all(found),'Input contains unmatched hemisphere/ROI names');
cb=colorbar(ax,'Location','southoutside');cb.Position=[.3 .035 .4 .025];publication_style(f);
parent=fileparts(output_stem);if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
savefig(f,[output_stem '.fig']);exportgraphics(f,[output_stem '.png'],'Resolution',180);close(f);
end

