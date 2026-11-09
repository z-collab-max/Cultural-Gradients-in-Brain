function reproduce_suppfig2(source_file,output_dir)
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));addpath(fullfile(root,'simulation'));addpath(fullfile(root,'figures','helpers'));
assert(nargin==2,'Supply source_file and output_dir');if ~isfolder(output_dir),mkdir(output_dir);end
d=load(source_file);
f=figure('Visible','off','Color','w','Position',[40 40 1500 1550]);
anchor=[.13 .40 .70 .94];rgb=[33 102 172;180 210 225;247 247 247;219 135 117;178 24 43]/255;
cm=interp1(linspace(0,1,5),rgb,linspace(0,1,256));
for j=1:2
 ax=axes(f,'Position',[.085 .765-(j-1)*.215 .82 .165]);a=squeeze(d.heatmaps(j,:,:));
 imagesc(ax,a,[-20 20]);ax.Color=[.85 .85 .85];h=findobj(ax,'Type','image');h.AlphaData=isfinite(a);colormap(ax,cm);
 hold(ax,'on');split=size(a,2)/2;xline(ax,split+.5,'k-');if size(a,1)>numel(d.common_centers),yline(ax,numel(d.common_centers)+.5,'k-');end
 set(ax,'YTick',1:numel(d.centers),'YTickLabel',string(d.centers),'XTick',[split/2 1.5*split],'XTickLabel',{'Left hemisphere','Right hemisphere'},'TickLength',[0 0],'TickLabelInterpreter','none');
 ylabel(ax,'Scanner center');title(ax,sprintf('Participant %d: spatially non-uniform deviations (%d centers)',j,d.scan_counts(j)),'FontSize',11,'FontWeight','normal');
 text(ax,1,1.08,sprintf('Median within-scan ROI span: %.1f percentage points',d.span_medians(j)),'Units','normalized','HorizontalAlignment','right','FontSize',8);
 cb=colorbar(ax);cb.Label.String='Thickness deviation (%)';
end
ax=axes(f,'Position',[.085 .29 .36 .19]);b=bar(ax,d.spatial_r,.65,'FaceColor','flat','EdgeColor','none');b.CData=repmat([.44 .38 .50],numel(d.spatial_r),1);b.CData(d.spatial_r<0,:)=repmat([.78 .47 .22],sum(d.spatial_r<0),1);
hold(ax,'on');yline(ax,0,'Color',[.5 .5 .5]);m=yline(ax,d.median_r,'--','Color',[.55 .55 .55]);m.Tag='MedianCorrelation';
set(ax,'XTick',1:numel(d.common_centers),'XTickLabel',string(d.common_centers),'XTickLabelRotation',65,'YLim',[-1 1],'TickDir','out');
ylabel(ax,'Cross-participant spatial correlation (r)');title(ax,'Regional patterns differ between participants','FontSize',11,'FontWeight','normal');
text(ax,.98,.92,sprintf('Median r = %.3f',d.median_r),'Units','normalized','HorizontalAlignment','right','FontSize',8);box(ax,'off');
ax=axes(f,'Position',[.66 .29 .27 .19]);hold(ax,'on');xline(ax,100,'--','Color',[.6 .6 .6]);colors=[33 102 172;180 91 80]/255;hh=gobjects(1,2);
for j=1:2
 y=(1:8)+(j-1.5)*.24;hh(j)=scatter(ax,d.relative_error(j,:),y,18,colors(j,:),'filled');
 for k=1:8,text(ax,d.relative_error(j,k)+.8,y(k),sprintf('%.1f%%',d.relative_error(j,k)),'Color',colors(j,:),'FontSize',7,'VerticalAlignment','middle');end
end
labels={'Raw','ROI-specific offset','Shared additive','Shared multiplicative','Shared linear','Shared quadratic','Shared cubic','Shared RBF'};
set(ax,'YDir','reverse','YTick',1:8,'YTickLabel',labels,'YLim',[.3 9.0],'XLim',[min(75,min(d.relative_error,[],'all')*.9) max(123,max(d.relative_error,[],'all')*1.05)],'TickDir','out');
xlabel(ax,{'Median absolute deviation remaining','relative to raw (= 100%)'});title(ax,'Cross-participant calibration at each center','FontSize',11,'FontWeight','normal');
legend(ax,hh,{'Participant 1','Participant 2'},'Location','southeast','Box','off');box(ax,'off');
ax=axes(f,'Position',[.05 .025 .9 .205]);axis(ax,[0 1 0 9]);axis(ax,'off');hold(ax,'on');
patch(ax,[0 1 1 0],[8 8 9 9],[.94 .94 .94],'EdgeColor','none');
text(ax,.01,8.5,'Method in panel c','FontSize',9,'FontWeight','bold');text(ax,.24,8.5,'Calibration formula and fitting','FontSize',9,'FontWeight','bold');
formula={'f(z) = z. No fit.','f_r(z_r) = z_r + a_r; a_r = y_r - x_r. ROI-specific offset.',...
 'f(z) = z + a; a = mean(y - x). Common offset.',...
 'f(z) = b z; b = sum(x y) / sum(x^2). Least squares through origin.',...
 'f(z) = a + b z. Ordinary least squares.',...
 'f(z) = a + b v + d v^2. Quadratic ridge fit.',...
 'f(z) = a + b v + d v^2 + e v^3. Cubic ridge fit.',...
 'f(z) = a + b v + sum_h w_h phi_h(v). Gaussian-basis ridge fit.'};
for j=1:8,text(ax,.01,8.5-j,labels{j},'FontSize',8,'Interpreter','none');text(ax,.24,8.5-j,formula{j},'FontSize',8,'Interpreter','none');end
plot(ax,[0 1],[9 9],'Color',[.5 .5 .5]);plot(ax,[0 1],[8 8],'Color',[.5 .5 .5]);plot(ax,[0 1],[0 0],'Color',[.5 .5 .5]);
positions=[.02 .94;.02 .49;.60 .49;.02 .235];
for j=1:4,annotation(f,'textbox',[positions(j,:) .03 .025],'String',char(96+j),'EdgeColor','none','FontSize',16);end
publication_style(f);savefig(f,fullfile(output_dir,'SuppFig2.fig'));exportgraphics(f,fullfile(output_dir,'SuppFig2.png'),'Resolution',180);close(f);
end
