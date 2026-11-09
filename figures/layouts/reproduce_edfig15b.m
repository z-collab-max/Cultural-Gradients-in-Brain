function reproduce_edfig15b(source_file,output_dir)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'simulation'));
assert(nargin==2,'Supply source_file and output_dir');if ~isfolder(output_dir),mkdir(output_dir);end
load(source_file);
set(groot,'defaultFigureVisible','off');
n_matched=size(score_joint_nations,1);
score_joint_provinces(:, 1) = -score_joint_provinces(:, 1);
score_joint_nations(:, 1)   = -score_joint_nations(:, 1);
corr_matrix(1, :)       = -corr_matrix(1, :);
corr_matrix_world(1, :) = -corr_matrix_world(1, :);
coeff_joint(:, 1)       = -coeff_joint(:, 1);
v_local_pc1 = [corr_matrix(1,1), corr_matrix(2,1)];
v_local_pc2 = [corr_matrix(1,2), corr_matrix(2,2)];
v_world_pc1 = [corr_matrix_world(1,1), corr_matrix_world(2,1)];
v_world_pc2 = [corr_matrix_world(1,2), corr_matrix_world(2,2)];
mu_china = mean(score_joint_provinces(:, 1:2), 1);
mu_world = mean(score_joint_nations(:, 1:2), 1);

scale = 2.4;
fig = figure('Position', [150, 150, 720, 480], 'Color', 'w');
ax = axes('Position', [0.12, 0.14, 0.54, 0.74]);
hold on;
plot_ellipse(score_joint_provinces(:, 1:2), '#E18F81');
plot_ellipse(score_joint_nations(:, 1:2), '#8FAADC');
s1 = scatter(score_joint_provinces(:,1), score_joint_provinces(:,2), 55, ...
    'MarkerFaceColor', '#E18F81', 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
s2 = scatter(score_joint_nations(:,1), score_joint_nations(:,2), 55, ...
    'MarkerFaceColor', '#8FAADC', 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
q1 = quiver(mu_china(1), mu_china(2), v_local_pc1(1)*scale, v_local_pc1(2)*scale, 0, ...
    'Color', '#D95319', 'LineWidth', 2.8, 'MaxHeadSize', 0.45);
q2 = quiver(mu_world(1), mu_world(2), v_world_pc1(1)*scale, v_world_pc1(2)*scale, 0, ...
    'Color', '#0072BD', 'LineWidth', 2.8, 'MaxHeadSize', 0.45);

q3 = quiver(mu_china(1), mu_china(2), v_local_pc2(1)*scale, v_local_pc2(2)*scale, 0, ...
    'Color', '#D95319', 'LineWidth', 2.8, 'MaxHeadSize', 0.45);
q4 = quiver(mu_world(1), mu_world(2), v_world_pc2(1)*scale, v_world_pc2(2)*scale, 0, ...
    'Color', '#0072BD', 'LineWidth', 2.8, 'MaxHeadSize', 0.45);
text(mu_china(1) + v_local_pc1(1)*scale + 0.15, mu_china(2) + v_local_pc1(2)*scale + 0.15, ...
    'Local PC1', 'Color', '#D95319', 'FontWeight', 'bold', 'FontSize', 9.5, 'HorizontalAlignment', 'left');
text(mu_world(1) + v_world_pc1(1)*scale - 0.15, mu_world(2) + v_world_pc1(2)*scale - 0.25, ...
    'World PC1', 'Color', '#0072BD', 'FontWeight', 'bold', 'FontSize', 9.5, 'HorizontalAlignment', 'center');
text(mu_china(1) + v_local_pc2(1)*scale - 0.2, mu_china(2) + v_local_pc2(2)*scale - 0.25, ...
    'Local PC2', 'Color', '#D95319', 'FontWeight', 'bold', 'FontSize', 9.5, 'HorizontalAlignment', 'right');
text(mu_world(1) + v_world_pc2(1)*scale + 0.15, mu_world(2) + v_world_pc2(2)*scale + 0.15, ...
    'World PC2', 'Color', '#0072BD', 'FontWeight', 'bold', 'FontSize', 9.5, 'HorizontalAlignment', 'left');
plot([-6 6], [0 0], 'k--', 'LineWidth', 0.5);
plot([0 0], [-6 6], 'k--', 'LineWidth', 0.5);
axis equal;
xlim([-5 5]);
ylim([-4 4]);
set(gca, 'XTick', -4:2:4, 'YTick', -4:2:4);

xlabel(sprintf('Joint PC1 (%.1f%%, flipped)', explained_joint(1)), 'FontSize', 10, 'FontWeight', 'bold');
ylabel(sprintf('Joint PC2 (%.1f%%)', explained_joint(2)), 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('PCA Subspace Geometry: China (%d provinces) vs. World (%d nations)',size(score_joint_provinces,1),n_matched), ...
    'FontSize', 10.5, 'FontWeight', 'bold');

lgd = legend([s1, s2, q1, q2], ...
    {sprintf('China Provinces (N=%d)',size(score_joint_provinces,1)), sprintf('World Nations (N=%d)', n_matched), 'China Local PC Axes', 'World PC Axes'}, ...
    'Position', [0.68, 0.62, 0.29, 0.25], 'Box', 'off', 'FontSize', 9);
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
ax.FontName = 'Helvetica';
ax.XLabel.FontName = 'Helvetica';
ax.YLabel.FontName = 'Helvetica';
ax.LineWidth = 0.85;
ax.XLabel.FontSize = 10;
ax.YLabel.FontSize = 10;
ax.TickLength = [0.02, 0.025];
box off;


out=output_dir;
publication_style(fig);savefig(fig,fullfile(out,'EDFig15b_joint_PCA_nations.fig'));
exportgraphics(fig,fullfile(out,'EDFig15b_joint_PCA_nations.png'),'Resolution',250);
print(fig,fullfile(out,'EDFig15b_joint_PCA_nations.svg'),'-dsvg','-painters');close(fig);
end
function plot_ellipse(data, color)
    mu = mean(data(:, 1:2), 1);
    Sigma = cov(data(:, 1:2));
    [V, D] = eig(Sigma);
    t = linspace(0, 2*pi, 100);
    circle = [cos(t); sin(t)];
    sd_scale = 2;
    ellipse = (V * sqrt(D) * sd_scale * circle)' + mu;
    plot(ellipse(:,1), ellipse(:,2), 'Color', color, 'LineWidth', 1.4, 'LineStyle', '-');
end

