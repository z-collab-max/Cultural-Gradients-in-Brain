function plot_matrix_results(input_file,output_file,kind)
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'simulation'));
D=load(input_file);f=figure('Visible','off','Color','w');
switch kind
 case 'scatter',scatter(D.x,D.y,20,'filled');xlabel('x');ylabel('y');
 case 'matrix',imagesc(D.matrix);axis image;colorbar;
 case 'summary',errorbar(D.x,D.mean,D.sem,'o');
 case 'surface'
  patch('Vertices',D.vertices,'Faces',D.faces,'FaceVertexCData',D.values,...
   'FaceColor','interp','EdgeColor','none');axis equal off;view(3);colorbar;
 otherwise,error('Unknown plot kind');
end
publication_style(f);folder=fileparts(output_file);if ~isempty(folder)&&~isfolder(folder),mkdir(folder);end
savefig(f,output_file);close(f);
end
