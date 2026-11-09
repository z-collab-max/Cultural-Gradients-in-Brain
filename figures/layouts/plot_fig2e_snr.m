function plot_fig2e_snr(source_file,output_dir)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'simulation'));
assert(nargin==2,'Supply source_file and output_dir');if ~isfolder(output_dir),mkdir(output_dir);end
d=load(source_file);
out=output_dir;
f=figure('Color','w','Position',[100 100 1200 510],'Visible','off');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for k=1:2
    ax=nexttile;hold(ax,'on');p=sprintf('pc%d_',k);
    x=d.([p 'x']);lo=d.([p 'lower']);hi=d.([p 'upper']);
    fill([x fliplr(x)],[lo fliplr(hi)],[.275 .49 .647],'FaceAlpha',.16,'EdgeColor','none','HandleVisibility','off');
    h1=plot(x,d.([p 'median']),'Color',[.275 .49 .647],'LineWidth',2);
    h2=plot(d.([p 'theory_x']),d.([p 'theory_y']),'--','Color',[.4 .4 .4],'LineWidth',1.4);
    q=d.([p 'actual_q']);
    h3=errorbar(d.([p 'actual_x']),q(2),q(2)-q(1),q(3)-q(2),'o','Color',[.71 .31 .28],'MarkerFaceColor',[.71 .31 .28],'CapSize',8);
    set(ax,'XScale','log','Box','off','XGrid','off','YGrid','off','FontSize',11);
    ylim([0 1.02]);tx=d.([p 'theory_x']);xlim([tx(1) tx(end)]);
    xlabel('Individual-level SNR (signal / noise variance)');title(sprintf('Cultural PC%d',k));
    if k==1
        ylabel('Recovery R^2');
        legend([h1 h2 h3],{'Calibrated mixed-model simulation','Moment approximation','Actual-SE-calibrated baseline'},'Location','southeast','Box','off','Color','none','FontSize',9);
    end
end
publication_style(f);savefig(f,fullfile(out,'individual_snr_recovery_editable.fig'));
exportgraphics(f,fullfile(out,'individual_snr_recovery_editable.png'),'Resolution',180);
close(f);
f=openfig(fullfile(out,'individual_snr_recovery_editable.fig'),'invisible');
assert(numel(findall(f,'Type','axes'))==2);close(f);
end
