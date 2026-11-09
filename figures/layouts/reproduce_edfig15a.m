function reproduce_edfig15a(source_file,output_dir)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'simulation'));
assert(nargin==2,'Supply source_file and output_dir');if ~isfolder(output_dir),mkdir(output_dir);end
load(source_file);
set(groot,'defaultFigureVisible','off');
vnames_short={'PDI','IDV','MAS','UAI','LTO','IVR'};display_names=string(world_nations_clean);prov_display=string(prov_display);n_prov=size(scores_china,1);n_world_clean=size(scores_world,1);
c_china    = '#E18F81';
c_asia     = '#9DC3E6';
c_eur_ame  = '#8FAADC';
c_loadings = '#D6604D';
fig = figure('Position', [509.8000  305.0000  720.0000  480.0000], 'Color', 'w');
ax = axes('Position', [0.08, 0.10, 0.72, 0.82]);
hold on;
plot([-8 8], [0 0], 'k--', 'LineWidth', 0.5);
plot([0 0], [-8 8], 'k--', 'LineWidth', 0.5);
plot_ellipse(scores_china, c_china, 2.0, 1.8);
idx_asia = find(continent_grp == 1);
plot_ellipse(scores_world(idx_asia, :), c_asia, 2.0, 1.8);
idx_eur = find(continent_grp == 2);
plot_ellipse(scores_world(idx_eur, :), c_eur_ame, 2.0, 1.8);
idx_ame = find(continent_grp == 3);
plot_ellipse(scores_world(idx_ame, :), [0.42, 0.55, 0.75], 2.0, 1.6);
h_china = scatter(scores_china(:,1), scores_china(:,2), 25, ...
    'MarkerFaceColor', c_china, 'MarkerEdgeColor', 'w', 'LineWidth', 0.1);
h_asia = scatter(scores_world(idx_asia, 1), scores_world(idx_asia, 2), 35, ...
    'MarkerFaceColor', c_asia, 'MarkerEdgeColor', [0.15 0.35 0.55], 'LineWidth', 0.1);
h_eur = scatter(scores_world(idx_eur, 1), scores_world(idx_eur, 2), 35, ...
    'MarkerFaceColor', c_eur_ame, 'MarkerEdgeColor', [0.20 0.30 0.50], 'LineWidth', 0.1);
h_ame = scatter(scores_world(idx_ame, 1), scores_world(idx_ame, 2), 35, ...
    'MarkerFaceColor', c_eur_ame, 'MarkerEdgeColor', [0.20 0.30 0.50], 'LineWidth', 0.1);
idx_oth = find(continent_grp == 4);
h_oth = scatter(scores_world(idx_oth, 1), scores_world(idx_oth, 2), 35, ...
    'MarkerFaceColor', c_eur_ame, 'MarkerEdgeColor', [0.20 0.30 0.50], 'LineWidth', 0.1);
scale_loadings = 3.6;
h_loading = [];
for v = 1:length(vnames_short)
    lx = coeff_joint(v, 1) * scale_loadings;
    ly = coeff_joint(v, 2) * scale_loadings;
    h_line = plot([0, lx], [0, ly], '--', 'Color', c_loadings, 'LineWidth', 1.3);
    if v == 1
        h_loading = h_line;
    end
    text(lx * 1.12, ly * 1.12, vnames_short{v}, 'Color', c_loadings, ...
        'FontSize', 9.5, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
fprintf('Placing non-overlapping country labels (no leader lines)...\n');
label_pos_world = compute_non_overlapping_positions_world(scores_world(:,1:2), display_names);

for i = 1:n_world_clean
    nm = display_names(i);
    lx = label_pos_world(i, 1);
    ly = label_pos_world(i, 2);
    
    if continent_grp(i) == 1
        t_col = [0.05, 0.30, 0.55];
    else
        t_col = [0.15, 0.25, 0.40];
    end
    
    text(lx, ly, nm, 'Color', t_col, 'FontSize', 6, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Helvetica');
end
fprintf('Placing non-overlapping English Chinese province labels...\n');
label_pos_china = compute_non_overlapping_positions_china(scores_china(:,1:2), prov_display);

for p = 1:n_prov
    pnm = prov_display(p);
    plx = label_pos_china(p, 1);
    ply = label_pos_china(p, 2);
    
    text(plx, ply, pnm, 'Color', [0.65, 0.18, 0.12], 'FontSize', 6,  ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Helvetica');
end
axis equal;
xlim([-6.0, 6]);
ylim([-4.5, 4.5]);
set(gca, 'XTick', -6:2:6, 'YTick', -4:2:4);

xlabel(sprintf('Joint PC1 (%.1f%% variance, flipped)', explained_joint(1)), ...
    'FontSize', 10.5, 'FontWeight', 'bold', 'FontName', 'Helvetica');
ylabel(sprintf('Joint PC2 (%.1f%% variance, flipped)', explained_joint(2)), ...
    'FontSize', 10.5, 'FontWeight', 'bold', 'FontName', 'Helvetica');
title(sprintf('Global Cultural Space: China (%d provinces) & World (%d nations)',n_prov,n_world_clean), ...
    'FontSize', 11.5, 'FontWeight', 'bold', 'FontName', 'Helvetica');
lgd = legend([h_china, h_asia, h_eur, h_ame, h_loading], ...
    {sprintf('China (N=%d)',n_prov),sprintf('Asia (N=%d)',numel(idx_asia)),sprintf('Europe (N=%d)',numel(idx_eur)),sprintf('Americas & Others (N=%d)',numel(idx_ame)+numel(idx_oth)),'Cultural dimensions'}, ...
    'Position', [0.79, 0.65, 0.20, 0.24], 'Box', 'off', 'FontSize', 8.5, 'FontName', 'Helvetica');
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
ax.FontName = 'Helvetica';
ax.LineWidth = 0.85;
ax.TickLength = [0.015, 0.02];
box off;


out=output_dir;
publication_style(fig);savefig(fig,fullfile(out,'EDFig15a_cultural_PCA.fig'));
exportgraphics(fig,fullfile(out,'EDFig15a_cultural_PCA.png'),'Resolution',250);
print(fig,fullfile(out,'EDFig15a_cultural_PCA.svg'),'-dsvg','-painters');close(fig);
end
function plot_ellipse(data, color, sd_scale, lw)
    if nargin < 4, lw = 1.5; end
    if size(data, 1) < 3, return; end
    mu = mean(data(:, 1:2), 1);
    Sigma = cov(data(:, 1:2));
    [V, D] = eig(Sigma);
    t = linspace(0, 2*pi, 100);
    circle = [cos(t); sin(t)];
    ellipse = (V * sqrt(D) * sd_scale * circle)' + mu;
    plot(ellipse(:,1), ellipse(:,2), 'Color', color, 'LineWidth', lw, 'LineStyle', '-');
end
function pos = compute_non_overlapping_positions_world(pts, names)
    n = size(pts, 1);
    pos = pts;
    for i = 1:n
        nm = names(i);
        switch nm
            case "Germany",        pos(i,:) = pts(i,:) + [-0.15, -0.28];
            case "Greece",         pos(i,:) = pts(i,:) + [ 0.42, -0.18];
            case "Italy",          pos(i,:) = pts(i,:) + [ 0.05,  0.30];
            case "Hungary",        pos(i,:) = pts(i,:) + [-0.05,  0.32];
            case "Czechia",        pos(i,:) = pts(i,:) + [ 0.42, -0.22];
            case "Slovakia",       pos(i,:) = pts(i,:) + [-0.35,  0.35];
            case "Canada",         pos(i,:) = pts(i,:) + [-0.48,  0.08];
            case "United States",  pos(i,:) = pts(i,:) + [-0.25,  0.35];
            case "Brazil",         pos(i,:) = pts(i,:) + [ 0.25,  0.32];
            case "Nigeria",        pos(i,:) = pts(i,:) + [ 0.42,  0.15];
            case "Malaysia",       pos(i,:) = pts(i,:) + [-0.38, -0.25];
            case "Singapore",      pos(i,:) = pts(i,:) + [-0.35, -0.25];
            case "Indonesia",      pos(i,:) = pts(i,:) + [ 0.52, -0.05];
            case "Japan",          pos(i,:) = pts(i,:) + [ 0.38,  0.08];
            case "Switzerland",    pos(i,:) = pts(i,:) + [-0.55, -0.18];
            case "Austria",        pos(i,:) = pts(i,:) + [-0.50,  0.25];
            case "South Africa",   pos(i,:) = pts(i,:) + [-0.55, -0.05];
            case "Colombia",       pos(i,:) = pts(i,:) + [-0.35, -0.28];
            case "Mexico",         pos(i,:) = pts(i,:) + [-0.35, -0.22];
            case "Uruguay",        pos(i,:) = pts(i,:) + [ 0.45, -0.18];
            case "Poland",         pos(i,:) = pts(i,:) + [-0.38, -0.15];
            case "Peru",           pos(i,:) = pts(i,:) + [ 0.35,  0.12];
            case "Romania",        pos(i,:) = pts(i,:) + [ 0.45,  0.08];
            case "Russia",         pos(i,:) = pts(i,:) + [ 0.38, -0.25];
            case "Latvia",         pos(i,:) = pts(i,:) + [ 0.35, -0.15];
            case "Ukraine",        pos(i,:) = pts(i,:) + [ 0.42, -0.12];
            case "Portugal",       pos(i,:) = pts(i,:) + [-0.42, -0.12];
            case "Spain",          pos(i,:) = pts(i,:) + [ 0.05, -0.28];
            case "France",         pos(i,:) = pts(i,:) + [ 0.42, -0.22];
            case "Netherlands",    pos(i,:) = pts(i,:) + [ 0.05, -0.25];
            case "Belgium",        pos(i,:) = pts(i,:) + [-0.05,  0.25];
            case "United Kingdom", pos(i,:) = pts(i,:) + [-0.35, -0.15];
            case "Norway",         pos(i,:) = pts(i,:) + [-0.38,  0.12];
            case "Sweden",         pos(i,:) = pts(i,:) + [-0.45, -0.12];
            case "Ireland",        pos(i,:) = pts(i,:) + [-0.42,  0.12];
            case "Denmark",        pos(i,:) = pts(i,:) + [-0.15, -0.28];
            case "Finland",        pos(i,:) = pts(i,:) + [-0.25, -0.25];
            case "Australia",      pos(i,:) = pts(i,:) + [-0.42, -0.18];
            case "New Zealand",    pos(i,:) = pts(i,:) + [-0.42, -0.22];
            case "India",          pos(i,:) = pts(i,:) + [-0.35,  0.20];
            case "Pakistan",       pos(i,:) = pts(i,:) + [-0.42,  0.15];
            case "Philippines",    pos(i,:) = pts(i,:) + [-0.15,  0.25];
            case "Thailand",       pos(i,:) = pts(i,:) + [ 0.45, -0.12];
            case "Vietnam",        pos(i,:) = pts(i,:) + [ 0.38,  0.15];
            case "South Korea",    pos(i,:) = pts(i,:) + [ 0.52,  0.12];
            case "Turkey",         pos(i,:) = pts(i,:) + [ 0.05,  0.28];
            otherwise,             pos(i,:) = pts(i,:) + [ 0.22,  0.22];
        end
    end
end
function pos = compute_non_overlapping_positions_china(pts, names)
    n = size(pts, 1);
    pos = pts;
    for i = 1:n
        nm = names(i);
        switch nm
            case "Inner Mongolia", pos(i,:) = pts(i,:) + [ 0.10,  0.22];
            case "Heilongjiang",   pos(i,:) = pts(i,:) + [ 0.30,  0.18];
            case "Jilin",          pos(i,:) = pts(i,:) + [ 0.20, -0.18];
            case "Liaoning",       pos(i,:) = pts(i,:) + [-0.25, -0.18];
            case "Beijing",        pos(i,:) = pts(i,:) + [-0.32,  0.22];
            case "Tianjin",        pos(i,:) = pts(i,:) + [-0.28, -0.15];
            case "Hebei",          pos(i,:) = pts(i,:) + [-0.22,  0.05];
            case "Shandong",       pos(i,:) = pts(i,:) + [-0.25, -0.20];
            case "Shanxi",         pos(i,:) = pts(i,:) + [-0.22,  0.18];
            case "Henan",          pos(i,:) = pts(i,:) + [ 0.25,  0.15];
            case "Jiangsu",        pos(i,:) = pts(i,:) + [ 0.28,  0.15];
            case "Zhejiang",       pos(i,:) = pts(i,:) + [ 0.25, -0.15];
            case "Anhui",          pos(i,:) = pts(i,:) + [ 0.25,  0.18];
            case "Jiangxi",        pos(i,:) = pts(i,:) + [ 0.25, -0.12];
            case "Fujian",         pos(i,:) = pts(i,:) + [ 0.22, -0.18];
            case "Guangdong",      pos(i,:) = pts(i,:) + [-0.28, -0.18];
            case "Guangxi",        pos(i,:) = pts(i,:) + [-0.25,  0.18];
            case "Hainan",         pos(i,:) = pts(i,:) + [ 0.25, -0.18];
            case "Hubei",          pos(i,:) = pts(i,:) + [-0.22,  0.18];
            case "Hunan",          pos(i,:) = pts(i,:) + [-0.25,  0.18];
            case "Sichuan",        pos(i,:) = pts(i,:) + [ 0.28,  0.15];
            case "Chongqing",      pos(i,:) = pts(i,:) + [-0.25, -0.18];
            case "Guizhou",        pos(i,:) = pts(i,:) + [-0.25,  0.18];
            case "Yunnan",         pos(i,:) = pts(i,:) + [ 0.25, -0.18];
            case "Shaanxi",        pos(i,:) = pts(i,:) + [ 0.22,  0.18];
            case "Gansu",          pos(i,:) = pts(i,:) + [-0.22, -0.15];
            case "Qinghai",        pos(i,:) = pts(i,:) + [ 0.28,  0.15];
            case "Ningxia",        pos(i,:) = pts(i,:) + [ 0.25, -0.18];
            case "Xinjiang",       pos(i,:) = pts(i,:) + [ 0.25,  0.18];
            case "Shanghai",       pos(i,:) = pts(i,:) + [-0.28, -0.18];
            otherwise,             pos(i,:) = pts(i,:) + [ 0.20,  0.15];
        end
    end
end
