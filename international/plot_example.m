function plot_example(output_dir)
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,output_dir=fullfile(root,'outputs','international');end
addpath(fullfile(root,'simulation'));
models={'PC1','PC2','USA_TL'};contexts={'endpoint'};
for m=1:numel(models)
 for c=1:numel(contexts)
  stem=[models{m} '_' contexts{c}];d=load(fullfile(output_dir,[stem '.mat']));
  f=figure('Visible','off','Color','w');scatter(d.targets(:,1),d.beta(:),24,'filled');
  xlabel('Synthetic cortical axis');ylabel('Cultural effect coefficient');
  publication_style(f);savefig(f,fullfile(output_dir,[stem '.fig']));close(f);
 end
end
end
