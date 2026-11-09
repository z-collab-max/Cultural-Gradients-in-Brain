function plot_geography(input_file,output_file)
% Input contains projected boundary polygons, xy, values and marker_sizes.
S=load(input_file);assert(size(S.xy,2)==2&&size(S.xy,1)==numel(S.values));
assert(all(isfinite(S.xy),'all')&&all(isfinite(S.values)));
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'simulation'));
f=figure('Color','w','Visible','off');ax=axes(f);hold(ax,'on');
for k=1:numel(S.boundaries)
 b=S.boundaries{k};patch(ax,b(:,1),b(:,2),[.95 .95 .95],'EdgeColor',[.65 .65 .65]);
end
scatter(ax,S.xy(:,1),S.xy(:,2),S.marker_sizes,S.values,'filled');axis(ax,'equal');axis(ax,'off');colorbar;
publication_style(f);parent=fileparts(output_file);if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
savefig(f,[output_file '.fig']);exportgraphics(f,[output_file '.png'],'Resolution',300);close(f);
end
