function reproduce_suppfig4(source_file,output_dir)
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));addpath(fullfile(root,'simulation'));addpath(fullfile(root,'figures','helpers'));
assert(nargin==2,'Supply source_file and output_dir');if ~isfolder(output_dir),mkdir(output_dir);end
d=load(source_file);
if iscell(d.route_labels),d.route_labels=cellfun(@char,d.route_labels,'UniformOutput',false);end
if ~isfolder(output_dir),mkdir(output_dir);end
f=figure('Visible','off','Color','w','Position',[50 50 1450 590]);
a=axes(f,'Position',[.065 .18 .36 .68]);bar(a,double(d.years),100*double(d.counts)/double(d.n),1,'FaceColor',[.53 .60 .71],'EdgeColor','w');hold(a,'on');
xline(a,d.median,'--','Color',[.35 .35 .35]);xlim(a,[-.5 max(d.years)+.5]);xticks(a,0:2:max(d.years));
xlabel(a,'Years since relocation');ylabel(a,'Participants (%)');title(a,sprintf('Relocation duration (n = %d)',d.n),'FontWeight','normal');
text(a,.97,.95,sprintf('Median = %g years\nIQR = %g–%g years',d.median,d.iqr(1),d.iqr(2)),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',10);
b=axes(f,'Position',[.64 .18 .32 .68]);bb=barh(b,1:numel(d.route_labels),d.route_percentages,'stacked','BarWidth',.72);
colors=[214 198 221;245 203 229;189 229 222;145 175 196]/255;
for k=1:4,bb(k).FaceColor=colors(k,:);bb(k).EdgeColor='none';end
set(b,'YDir','reverse','YTick',1:numel(d.route_labels),'YTickLabel',string(d.route_labels),'TickLabelInterpreter','none');
totals=sum(d.route_percentages,2);xlim(b,[0 max(totals)*1.18]);ylim(b,[.3 numel(d.route_labels)+.7]);
for k=1:numel(d.route_labels),text(b,totals(k)+max(totals)*.015,k,sprintf('%.2f%%',totals(k)),'FontSize',8,'VerticalAlignment','middle');end
xlabel(b,sprintf('Participants (%% of all %d participants)',d.n));title(b,'Most frequent relocation routes','FontWeight','normal');
legend(b,bb,{'0–1 years','2–3 years','4–5 years','≥6 years'},'Location','southeast','NumColumns',2,'Box','off');
for ax=[a b],set(ax,'Box','off','TickDir','out');end
annotation(f,'textbox',[.015 .89 .03 .05],'String','a','EdgeColor','none','FontSize',16);
annotation(f,'textbox',[.49 .89 .03 .05],'String','b','EdgeColor','none','FontSize',16);
publication_style(f);savefig(f,fullfile(output_dir,'SuppFig4.fig'));exportgraphics(f,fullfile(output_dir,'SuppFig4.png'),'Resolution',180);close(f);
end
