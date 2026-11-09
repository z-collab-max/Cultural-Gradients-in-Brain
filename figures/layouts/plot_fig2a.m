function plot_fig2a(source_file,out)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'simulation'));addpath(fullfile(root,'figures','helpers'));
assert(nargin==2);if ~isfolder(out),mkdir(out);end
set(groot,'defaultFigureVisible','off');load(source_file);
d = Fig2a_Data;
d.score(:,2)=-d.score(:,2);
d.coeff(:,2)=-d.coeff(:,2);
figure('Color', 'white', 'Position', [100 200 900 600], 'Name', 'PCA Biplot');
scatter(d.score(:,1), d.score(:,2), 30, 'filled', ...
    'MarkerFaceColor', hex2rgb('#71819E'), 'MarkerEdgeColor', 'none'); 
hold on;
text(d.score(:,1) + 0.1, d.score(:,2), d.prov_use, ...
    'FontSize', 5, 'Color', [0.5 0.5 0.5], 'Interpreter', 'none');
max_score = max(abs(d.score(:)));
max_coeff = max(abs(d.coeff(:)));
scaleFactor = (max_score / max_coeff) * 0.7; 

for i = 1:size(d.coeff, 1)
    x_v = d.coeff(i,1) * scaleFactor;
    y_v = d.coeff(i,2) * scaleFactor;
    line([0 x_v], [0 y_v], 'Color', hex2rgb('#B88076'), 'LineWidth', 0.85); 
    text(x_v * 1.1, y_v * 1.1, d.var_labels_short{i}, ...
        'Color', hex2rgb('#AB4D3E'), 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end
xlabel(sprintf('PC1 (%.2f%%)', d.explained(1)), 'FontSize', 8);
ylabel(sprintf('PC2 (%.2f%%)', d.explained(2)), 'FontSize', 8);
title('\fontsize{11}\rm PCA Biplot: provinces and culture loadings');

axis equal; grid off;
xlim([-4 7]);ylim([-4 3]);
set(gca, 'TickDir', 'out', 'LineWidth', 0.85, 'FontSize', 8, 'TickLength', 4*[0.01 0.025]);
x_lim = xlim; y_lim = ylim;
plot(x_lim, [0 0], 'k--', 'LineWidth', 0.85, 'Color', [0.4 0.4 0.4]);
plot([0 0], y_lim, 'k--', 'LineWidth', 0.85, 'Color', [0.4 0.4 0.4]);
hold off;
figure('Color', 'white', 'Position', 0.8*[100 200 550 300], 'Name', 'Scree Plot');
yyaxis left
bar_handle = bar(d.explained);
bar_handle.FaceColor = hex2rgb('#B3C1D9');
bar_handle.EdgeColor = 'none';
ylabel('Variance (%)');
xlabel('PCs');
set(gca, 'YColor', 'k');
text(1:length(d.explained), d.explained, num2str(d.explained, '%.1f%%'), ...
    'vert', 'bottom', 'horiz', 'center', 'FontSize', 8);
yyaxis right
plot(cumsum(d.explained), '-o', 'LineWidth', 0.85, 'Color', hex2rgb('#B88076'), ...
    'MarkerFaceColor', 'w', 'MarkerSize', 4);
ylabel('Cumulative Variance (%)');
ylim([0 110]);
set(gca, 'YColor', hex2rgb('#AB4D3E'));
title('\fontsize{11}\rm Variance explained across provinces');
box off;
set(gca, 'TickDir', 'out', 'LineWidth', 0.85, 'FontSize', 8, 'TickLength', 4*[0.01 0.025]);


figs=findall(groot,'Type','figure');
for ii=1:numel(figs)
 publication_style(figs(ii));exportgraphics(figs(ii),fullfile(out,sprintf('plot_fig2a_%d.png',ii)),'Resolution',180);
 publication_style(figs(ii));savefig(figs(ii),fullfile(out,sprintf('plot_fig2a_%d.fig',ii)));
end
close all;
end
