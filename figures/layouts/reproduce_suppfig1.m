function reproduce_suppfig1(source_file,output_dir)
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));addpath(fullfile(root,'simulation'));addpath(fullfile(root,'figures','helpers'));
assert(nargin==2,'Supply source_file and output_dir');if ~isfolder(output_dir),mkdir(output_dir);end
d=load(source_file);
d.dictionary=cellfun(@char,d.dictionary,'UniformOutput',false);
d.word_labels=cellfun(@char,d.word_labels,'UniformOutput',false);
f=figure('Visible','off','Color','w','Position',[50 50 1500 1000]);
ax=axes(f,'Position',[.025 .57 .34 .36]);axis(ax,[0 4 0 11]);axis(ax,'off');hold(ax,'on');
patch(ax,[0 4 4 0],[10 10 11 11],[.12 .27 .46],'EdgeColor','none');
text(ax,1,10.5,'tightness','Color','w','HorizontalAlignment','center','FontSize',11);
text(ax,3,10.5,'looseness','Color','w','HorizontalAlignment','center','FontSize',11);
for row=1:10
 y=10-row;
 if mod(row,2)==0,patch(ax,[0 4 4 0],[y y y+1 y+1],[.95 .96 .98],'EdgeColor','none');end
 for col=1:4
  h=text(ax,col-.5,y+.5,d.dictionary{row,col},'HorizontalAlignment','center','Interpreter','none','FontSize',9);
  if ismember(col,[2 4]),h.Tag='ChineseWord';end
 end
end
plot(ax,[2 2],[0 10],'Color',[.7 .7 .7]);
for j=1:2
 ax=axes(f,'Position',[.43+(j-1)*.29 .59 .20 .32]);e=d.explained(j,:);
 yyaxis(ax,'left');bar(ax,1:10,e,1,'FaceColor',[.85 .89 .95],'EdgeColor','none');ylabel(ax,'Variance (%)');
 yyaxis(ax,'right');plot(ax,1:10,cumsum(e),'-o','Color',[.72 .50 .46],'MarkerSize',3,'MarkerFaceColor',[.72 .50 .46]);ylabel(ax,'Cumulative variance (%)');ylim(ax,[0 105]);
 ax.YAxis(1).Color=[.2 .2 .2];ax.YAxis(2).Color=[.67 .30 .24];xticks(ax,1:10);xlim(ax,[.3 10.7]);xlabel(ax,'Principal components');box(ax,'off');
 titles={'Tightness PCA','Looseness PCA'};title(ax,titles{j},'FontWeight','normal','FontSize',11);
end
colors=[.18 .45 .71;.85 .33 .10;.12 .63 .53];
titles={'Composite difference (PC_T1 - PC_L1)','Pure tightness (PC_T1)','Pure looseness (PC_L1)'};
labels={'Weibo tightness difference (normalized)','Weibo tightness component (normalized)','Weibo looseness component (normalized)'};
for j=1:3
 ax=axes(f,'Position',[.07+(j-1)*.32 .09 .25 .34]);hold(ax,'on');a=squeeze(d.fit(j,:,:));c=colors(j,:);
 fill(ax,[a(:,1);flipud(a(:,1))],[a(:,3);flipud(a(:,4))],c,'FaceAlpha',.18,'EdgeColor','none');
 plot(ax,a(:,1),a(:,2),'Color',c);scatter(ax,d.x(:,j),d.y,24,c,'filled','MarkerEdgeColor','w');
 for n=1:numel(d.y)
  dx=.02*range(d.x(:,j));align='left';if d.x(n,j)>1.8,dx=-dx;align='right';end
  text(ax,d.x(n,j)+dx,d.y(n),d.province_names{n},'HorizontalAlignment',align,'FontSize',6,'Color',[.3 .3 .3],'Interpreter','none');
 end
 text(ax,.05,.92,sprintf('r = %.3f, p = %.3f',d.r(j),d.p(j)),'Units','normalized','FontSize',9);
 title(ax,titles{j},'Interpreter','tex','FontSize',10,'FontWeight','normal');xlabel(ax,labels{j});ylabel(ax,'TL score');ylim(ax,[-.25 .25]);xlim(ax,[-2.4 2.4]);box(ax,'off');set(ax,'TickDir','out');
end
for item={{'a',[.015 .94 .03 .04]},{'b',[.39 .94 .03 .04]},{'c',[.015 .45 .03 .04]}}
 q=item{1};annotation(f,'textbox',q{2},'String',q{1},'EdgeColor','none','FontSize',16);
end
if ~isfolder(output_dir),mkdir(output_dir);end
publication_style(f);
for h=reshape(findall(f,'Tag','ChineseWord'),1,[]),h.FontName='Microsoft YaHei';end
drawnow;savefig(f,fullfile(output_dir,'SuppFig1.fig'));
print(f,fullfile(output_dir,'SuppFig1.svg'),'-dsvg','-painters');
exportgraphics(f,fullfile(output_dir,'SuppFig1.png'),'Resolution',200);close(f);
end
