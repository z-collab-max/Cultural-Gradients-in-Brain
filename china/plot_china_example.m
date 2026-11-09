function plot_china_example(output_dir)
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'simulation'));
if nargin<1,output_dir=fullfile(root,'outputs','china');end
source=fullfile(output_dir,'source_data');
f=figure('Visible','off','Color','w');tiledlayout(1,2);
for pc=1:2
 d=readtable(fullfile(source,sprintf('mixed_PC%d.csv',pc)));nexttile;
 errorbar(d.roi,d.beta,d.se,'o');xlabel('ROI');ylabel('Standardized coefficient');title(sprintf('PC%d',pc));
end
publication_style(f);savefig(f,fullfile(output_dir,'cultural_effects.fig'));close(f);
f=figure('Visible','off','Color','w');tiledlayout(1,2);
for pc=1:2
 d=readtable(fullfile(source,sprintf('relocation_PC%d.csv',pc)));nexttile;
 plot(d.roi,[d.stayer,d.relocator],'-o');xlabel('ROI');ylabel('Standardized coefficient');
 title(sprintf('PC%d',pc));legend('Stayer','Relocator','Location','best');
end
publication_style(f);savefig(f,fullfile(output_dir,'relocation_effects.fig'));close(f);
files={'behavior_RDM.mat','mediation_RDM_PC1.mat','mediation_RDM_PC2.mat'};
for j=1:numel(files)
 d=load(fullfile(source,files{j}));o=d.order;f=figure('Visible','off','Color','w');tiledlayout(1,2);
 nexttile;imagesc(o.model_rdm(o.best_order,o.best_order));axis square;colorbar;
 xticks(1:7);yticks(1:7);xticklabels(o.task_names);yticklabels(o.task_names);xtickangle(45);title('Selected task RDM');
 nexttile;bar(o.rsa_map);xlabel('ROI');ylabel('Spearman RSA coefficient');
 publication_style(f);[~,stem]=fileparts(files{j});savefig(f,fullfile(output_dir,[stem '.fig']));close(f);
end
end
