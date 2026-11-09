function plot_fig2e_replacement(source_file,out)
here=fileparts(mfilename('fullpath'));addpath(fullfile(fileparts(fileparts(here)),'simulation'));
if ~isfolder(out),mkdir(out);end
d=load(source_file);d=replacement_significance(d);
red=[.71 .31 .28];grey=[.47 .47 .47];
f=figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1100 470]);
for k=1:2
    ax=subplot(1,2,k);hold(ax,'on');p=sprintf('pc%d_',k);x=d.([p 'x']);
    tl=d.([p 'top_lower']);tm=d.([p 'top_median']);tu=d.([p 'top_upper']);rl=d.([p 'random_lower']);rm=d.([p 'random_median']);ru=d.([p 'random_upper']);
    fill(ax,[x fliplr(x)],[tl fliplr(tu)],red,'FaceAlpha',.14,'EdgeColor','none','HandleVisibility','off');
    fill(ax,[x fliplr(x)],[rl fliplr(ru)],grey,'FaceAlpha',.14,'EdgeColor','none','HandleVisibility','off');
    h1=plot(ax,x,tm,'Color',red,'LineWidth',2);h2=plot(ax,x,rm,'Color',grey,'LineWidth',2);sig=d.([p 'significant_pct']);
    if ~isempty(sig)
        line(ax,[min(sig) max(sig)],[.975 .975],'Color',[.14 .14 .14],'LineWidth',1.4,'HandleVisibility','off');
        text(ax,mean([min(sig) max(sig)]),1.005,'Holm p < 0.05','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    end
    title(ax,sprintf('Cultural PC%d',k));xlabel(ax,'ROIs replaced by signed median (%)');
    if k==1
        ylabel(ax,'Recovery r: estimate vs modified target beta map');lg=legend(ax,[h1 h2],{'Largest observed |t| first','Random same-number ROIs'},'Location','southwest');set(lg,'Box','off','Color','none');
    end
    xlim(ax,[0 99]);ylim(ax,[-.2 1.06]);set(ax,'Box','off','FontName','Arial','FontSize',11,'TickDir','out');grid(ax,'off');
end
publication_style(f);savefig(f,fullfile(out,'province_actual_se_roi_ablation_editable.fig'));
exportgraphics(f,fullfile(out,'province_actual_se_roi_ablation_editable.png'),'Resolution',180);
close(f);
end

function d=replacement_significance(d)
p=[];for k=1:2,pre=sprintf('pc%d_',k);v=d.([pre 'null']);b=d.([pre 'observed']);p=[p mean(v<=b(:)',1)];end
[sorted,index]=sort(p);adjusted=zeros(size(p));running=0;
for j=1:numel(p),running=max(running,min(1,(numel(p)-j+1)*sorted(j)));adjusted(index(j))=running;end
for k=1:2,pre=sprintf('pc%d_',k);pct=d.([pre 'test_pct']);n=numel(pct);sig=adjusted((k-1)*n+(1:n))<.05;d.([pre 'significant_pct'])=pct(sig);end
end
