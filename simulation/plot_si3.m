function plot_si3(data_dir)
here=fileparts(mfilename('fullpath'));P=fileparts(here);addpath(here);
if nargin<1,data_dir=fullfile(P,'simulation','reference','si3');end
f=figure('Visible','off','Color','w','Position',[40 40 2700 1250]);
colors=[107 152 187;105 170 162;151 143 184;209 168 102;217 135 124]/255;ref=[184 189 193]/255;
px=[.04 .15 .28 .39 .54 .65 .78 .89];w=.075;h=w*2700/1250;ks=[1 3 5];
for row=1:3
 y=.72-(row-1)*.30;
 labelax=axes(f,'Position',[.001 y .012 h],'Visible','off');
 text(labelax,.5,.5,sprintf('%d artifact sources',ks(row)),'Rotation',90,'HorizontalAlignment','center','FontSize',10);
 for arm=1:2
  arms={'M','C'};d=load(fullfile(data_dir,sprintf('SuppFig3_K%d_%s.mat',ks(row),arms{arm})));
  for col=1:4
   ax=axes(f,'Position',[px((arm-1)*4+col) y w h]);hold(ax,'on');
   if col<=2
    x=d.example_gt;if col==1,v=d.example_raw;else,v=d.example_corrected;end;limits=[-4 4];
    xlabel(ax,'ground truth coefficient (z-score)');if col==1,ylabel(ax,{'estimated coefficient','(z-score)'});end
   else
    x=d.true_gradient_r;if col==3,v=d.raw_gradient_r;else,v=d.corrected_gradient_r;end;limits=[-1.1 1.1];
    xlabel(ax,'r (cultural gradient, SA)');if col==3,ylabel(ax,'r (cultural gradient, SA)');end
    plot(ax,limits,[0 0],':','Color',ref);plot(ax,[0 0],limits,':','Color',ref);
   end
   plot(ax,limits,limits,'--','Color',ref);
   xx=[min(x(:)) max(x(:))];mid=mean(xx);xx=mid+[-1 1]*diff(xx)*1.25/2;
   for b=1:5,p=polyfit(x(b,:),v(b,:),1);plot(ax,xx,polyval(p,xx),'Color',colors(b,:));end
   hh=gobjects(5,1);
   for b=1:5,hh(b)=scatter(ax,x(b,:),v(b,:),7,colors(b,:),'filled','MarkerFaceAlpha',.60,'MarkerEdgeColor','none');end
   if row==1&&arm==1&&col==3,lg=legend(ax,hh,string(d.artifact_rho),'Location','northwest','Box','off');lg.Title.String='r (artifact, SA)';end
   if mod(col,2)==1,title(ax,'Raw','FontSize',10,'FontWeight','normal');else,title(ax,'De-artifacted','FontSize',10,'FontWeight','normal');end
   xlim(ax,limits);ylim(ax,limits);axis(ax,'square');box(ax,'off');set(ax,'TickDir','out');
  end
 end
end
for arm=1:2
 x=.025+(arm-1)*.5;labels={'a','b'};heading={'Overlap = 80%, threshold = 100%','Overlap = 100%, threshold = 99.9%'};
 annotation(f,'textbox',[x .95 .03 .04],'String',labels{arm},'EdgeColor','none','FontSize',18);
 annotation(f,'textbox',[x+.04 .95 .41 .04],'String',heading{arm},'EdgeColor','none','HorizontalAlignment','center','FontSize',12);
 annotation(f,'textbox',[x+.02 .90 .20 .04],'String','Cultural gradient recovery','EdgeColor','none','HorizontalAlignment','center','FontSize',11);
 annotation(f,'textbox',[x+.24 .90 .23 .04],'String','r (cultural gradient, S-A axis) recovery','EdgeColor','none','HorizontalAlignment','center','FontSize',11);
end
save_supplementary_figure(f,'SupplementaryFigure3','SuppFig3');close(f);
end
