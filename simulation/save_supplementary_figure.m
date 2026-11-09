function save_supplementary_figure(f,folder,stem)
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'outputs','simulation_figures',folder);if ~isfolder(out),mkdir(out);end
publication_style(f);
for h=reshape(findall(f,'Tag','ChineseWord'),1,[]),h.FontName='Microsoft YaHei';end
drawnow;savefig(f,fullfile(out,[stem '.fig']));
print(f,fullfile(out,[stem '.svg']),'-dsvg','-painters');
exportgraphics(f,fullfile(out,[stem '.png']),'Resolution',200);
fprintf('EXPORTED %s/%s\n',folder,stem);
end
