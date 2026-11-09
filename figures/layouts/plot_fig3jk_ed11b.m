function plot_fig3jk_ed11b(source_file,out)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'simulation'));addpath(fullfile(root,'figures','helpers'));
assert(nargin==2);if ~isfolder(out),mkdir(out);end
previous=get(groot,'defaultFigureVisible');cleanup=onCleanup(@()set(groot,'defaultFigureVisible',previous));existing=findall(groot,'Type','figure');
set(groot,'defaultFigureVisible','off');load(source_file);
for ccc = 1:2
    d = Fig4_ED9_Data(ccc);
    n_pos = length(d.idx_pos_sorted);
    n_neg = length(d.idx_neg_sorted);
    d.y_pos_list = (n_pos + n_neg : -1 : n_neg + 1)';
    d.y_neg_list = (n_neg : -1 : 1)';

    pc_title = sprintf('Culture PC%d', ccc);
    figure('Name', sprintf('%s - Effect Size Dist', pc_title), 'Color', 'w');
    set(gcf,'Position',[100 400 0.8*0.8*450 0.8*0.8*200]);
    hold on;
    histogram(d.cohendpool_control_posM, 'Normalization', 'probability', 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'FaceColor', colorpair{1,2});
    histogram(d.cohendpool_control_negM, 'Normalization', 'probability', 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'FaceColor', colorpair{2,2});
    line([d.treatcohendMN d.treatcohendMN], [0 0.06], 'LineWidth', 0.85, 'Color', colorpair{2,1});
    line([d.treatcohendMP d.treatcohendMP], [0 0.06], 'LineWidth', 0.85, 'Color', colorpair{1,1});
    line([d.treatcohend_MHH_MLL_N d.treatcohend_MHH_MLL_N], [0 0.06], 'LineWidth', 0.85, 'Color', colorpair{2,3}, 'LineStyle', '--');
    line([d.treatcohend_MHH_MLL_P d.treatcohend_MHH_MLL_P], [0 0.06], 'LineWidth', 0.85, 'Color', colorpair{1,3}, 'LineStyle', '--');

    title(pc_title, 'FontSize', 10);
    xlabel("Cohen's d (High - Low)", 'FontSize', 8); ylabel('Probability', 'FontSize', 8);
    set(gca, 'XLim', [-0.5 0.5], 'TickDir', 'out', 'Box', 'off', 'FontSize', 8, 'LineWidth', 0.85);
    figure('Name', sprintf('%s - Interaction', pc_title), 'Color', 'w');
    set(gcf,'Position',[200 300 280 220]);
    hold on;
    histogram(d.diff_dist, 'Normalization', 'pdf', 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
    mu = mean(d.diff_dist); sigma = std(d.diff_dist);
    x_lims = get(gca, 'XLim');
    plot_min = min([x_lims(1), d.diff_treat_cross, d.diff_treat_same]) * 1.1;
    plot_max = max([x_lims(2), d.diff_treat_cross, d.diff_treat_same]) * 1.1;
    x_fit = linspace(plot_min, plot_max, 100);
    plot(x_fit, normpdf(x_fit, mu, sigma), 'Color', [0.3 0.3 0.3], 'LineWidth', 0.85);
    ylim_curr = get(gca, 'YLim'); y_max_line = ylim_curr(2);
    line([d.diff_treat_cross d.diff_treat_cross], [0 y_max_line], 'Color', '#b2182b', 'LineWidth', 0.85, 'LineStyle', '-');
    line([d.diff_treat_same d.diff_treat_same], [0 y_max_line], 'Color', '#d6604d', 'LineWidth', 0.85, 'LineStyle', '--');
    title({sprintf('\\fontsize{8}\\rm%s', pc_title), ...
        ['\fontsize{6}\rm Cross p=' sprintf('%.3f', d.p_val_cross) ', Same p=' sprintf('%.3f', d.p_val_same)]}, 'Interpreter', 'tex');
    xlabel(['\fontsize{5} \rm Difference of (Pos-Neg) Differences' newline '[(Pos-Neg)_{High} - (Pos-Neg)_{Low}]']);
    xlim([plot_min, plot_max]);
    set(gca, 'TickDir', 'out', 'LineWidth', 0.85, 'FontSize', 8, 'TickLength', 4*[0.01 0.025], 'YColor', 'none', 'Box', 'off');
    figure('Name', sprintf('%s - Single Axis Forest', pc_title), 'Color', 'w', 'Position', [100 400 0.8*0.8*450 0.9*400]);
    hold on;

    c_pos_stable = colorpair{1,3}; c_pos_trans = colorpair{1,2}; c_pos_line = '#CAA19A';
    c_neg_stable = colorpair{2,3}; c_neg_trans = colorpair{2,2}; c_neg_line = '#99A4B9';

    line_data_pos_x = d.mean_stable_pos(d.idx_pos_sorted); line_data_pos_y = d.y_pos_list(:); line_data_pos_trans = d.obs_d_pos(d.idx_pos_sorted);
    line_data_neg_x = d.mean_stable_neg(d.idx_neg_sorted); line_data_neg_y = d.y_neg_list(:); line_data_neg_trans = d.obs_d_neg(d.idx_neg_sorted);
    fill([line_data_pos_x; flipud(line_data_pos_trans)], [line_data_pos_y; flipud(line_data_pos_y)], hex2rgb(c_pos_trans), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([line_data_neg_x; flipud(line_data_neg_trans)], [line_data_neg_y; flipud(line_data_neg_y)], hex2rgb(c_neg_trans), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(line_data_pos_x, line_data_pos_y, '-', 'Color', hex2rgb(c_pos_stable), 'LineWidth', 0.85);
    plot(line_data_neg_x, line_data_neg_y, '-', 'Color', hex2rgb(c_neg_stable), 'LineWidth', 0.85);

    cap_size = 0.25; lw_ci = 1;
    for i = 1:length(d.idx_pos_sorted)
        y = d.y_pos_list(i); curr_idx = d.idx_pos_sorted(i);
        ci = d.ci_stable_pos(curr_idx, :); val_stable = d.mean_stable_pos(curr_idx); val_trans = d.obs_d_pos(curr_idx);

        col = hex2rgb(c_pos_line);
        plot(ci, [y y], '-', 'Color', col, 'LineWidth', lw_ci);
        plot([ci(1) ci(1)], [y-cap_size, y+cap_size], '-', 'Color', col, 'LineWidth', lw_ci);
        plot([ci(2) ci(2)], [y-cap_size, y+cap_size], '-', 'Color', col, 'LineWidth', lw_ci);

        scatter(val_stable, y, 30, 'filled', 'MarkerFaceColor', c_pos_stable, 'MarkerEdgeColor', 'none');
        scatter(val_trans, y, 30, 'filled', 'MarkerFaceColor', c_pos_trans, 'MarkerEdgeColor', 'none');
    end
    for i = 1:length(d.idx_neg_sorted)
        y = d.y_neg_list(i); curr_idx = d.idx_neg_sorted(i);
        ci = d.ci_stable_neg(curr_idx, :); val_stable = d.mean_stable_neg(curr_idx); val_trans = d.obs_d_neg(curr_idx);

        col = hex2rgb(c_neg_line);
        plot(ci, [y y], '-', 'Color', col, 'LineWidth', lw_ci);
        plot([ci(1) ci(1)], [y-cap_size, y+cap_size], '-', 'Color', col, 'LineWidth', lw_ci);
        plot([ci(2) ci(2)], [y-cap_size, y+cap_size], '-', 'Color', col, 'LineWidth', lw_ci);

        scatter(val_stable, y, 30, 'filled', 'MarkerFaceColor', c_neg_stable, 'MarkerEdgeColor', 'none');
        scatter(val_trans, y, 30, 'filled', 'MarkerFaceColor', c_neg_trans, 'MarkerEdgeColor', 'none');
    end
    ylim([0, max(d.y_pos_list)+2]); xline(0, '--k', 'Alpha', 0.5);
    text(0, max(d.y_pos_list)+1, sprintf('\\leftarrow Positive ROIs\n'), 'Color', c_pos_stable, 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    text(0, max(d.y_neg_list)+1, sprintf('Negative ROIs \\rightarrow\n '), 'Color', c_neg_stable, 'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    title(pc_title, 'FontSize', 10);
    xlabel("Cohen's d (Effect Size)", 'FontSize', 8);
    set(gca, 'YColor', 'none', 'XLim', [-0.4 0.4], 'XTick', -0.4:0.2:0.4, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 0.85, 'FontSize', 8, 'TickLength', 4*[0.01 0.025]);
    figure('Name', sprintf('%s - Dual Panel Forest', pc_title), 'Color', 'w', 'Position', [100, 100, 400, 400]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    c_LL_stable = colorpair{2,3}; c_LH_trans = colorpair{2,2}; c_Low_area = colorpair{2,1};
    c_HH_stable = colorpair{1,3}; c_HL_trans = colorpair{1,2}; c_High_area = colorpair{1,1};
    ax1 = nexttile; hold on;
    title({sprintf('\\fontsize{8}\\rm%s', pc_title), '\fontsize{11}\rm Low Origin (LL \rightarrow LH)'}, 'FontSize', 11, 'Color', hex2rgb(c_LL_stable));
    plot_forest_part(d.y_pos_list, d.dp_data_pos(:,1:4), c_LL_stable, c_LH_trans, c_Low_area);
    plot_forest_part(d.y_neg_list, d.dp_data_neg(:,1:4), c_LL_stable, c_LH_trans, c_Low_area);
    xline(0, '--k'); xlabel('\fontsize{9}\rm Z-scored Thickness');
    ylim([0, max(d.y_pos_list)+2]); xlim([-0.2 0.2]);
    set(gca, 'YColor', 'none', 'TickDir', 'out', 'Box', 'off', 'FontSize', 8);
    text(0, max(d.y_pos_list)+1, '\fontsize{6}\rm Positive ROIs', 'Horiz','center', 'Color','k');
    text(0, max(d.y_neg_list)+1, '\fontsize{6}\rm Negative ROIs', 'Horiz','center', 'Color','k');
    ax2 = nexttile; hold on;
    title({sprintf('\\fontsize{8}\\rm%s', pc_title), '\fontsize{11}\rm High Origin (HH \rightarrow HL)'}, 'Color', hex2rgb(c_HH_stable));
    plot_forest_part(d.y_pos_list, d.dp_data_pos(:,5:8), c_HH_stable, c_HL_trans, c_High_area);
    plot_forest_part(d.y_neg_list, d.dp_data_neg(:,5:8), c_HH_stable, c_HL_trans, c_High_area);
    xline(0, '--k'); xlabel('\fontsize{9}\rm Z-scored Thickness');
    ylim([0, max(d.y_pos_list)+2]); xlim([-0.2 0.2]);
    set(gca, 'YColor', 'none', 'TickDir', 'out', 'Box', 'off', 'FontSize', 8);
    text(0, max(d.y_pos_list)+1, '\fontsize{6}\rm Positive ROIs', 'Horiz','center', 'Color','k');
    text(0, max(d.y_neg_list)+1, '\fontsize{6}\rm Negative ROIs', 'Horiz','center', 'Color','k');
    linkaxes([ax1, ax2], 'xy');
end
if ~exist(out,'dir'),mkdir(out);end
figs=setdiff(findall(groot,'Type','figure'),existing,'stable');
for j=1:numel(figs)
 set(figs(j),'Position',[100 100 1050 650]);
 name=regexprep(get(figs(j),'Name'),'[^a-zA-Z0-9_]','_');if isempty(name),name=sprintf('panel%d',j);end
 publication_style(figs(j));exportgraphics(figs(j),fullfile(out,['plot_fig3jk_ed11b_' name '.png']),'Resolution',150);
 savefig(figs(j),fullfile(out,['plot_fig3jk_ed11b_' name '.fig']));
end
close(figs);disp('plot_fig3jk_ed11b DONE');
end

function plot_forest_part(y_list, data, c_stable, c_move, c_area)
y_list = y_list(:);
x_stable = data(:,1); x_move = data(:,2);
ci_L = data(:,3); ci_H = data(:,4);

X_poly = [x_stable; flipud(x_move)]; Y_poly = [y_list; flipud(y_list)];
fill(X_poly, Y_poly, hex2rgb(c_area), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

cap = 0.3; lw = 1;
for i = 1:length(y_list)
    yy = y_list(i); line_col = [hex2rgb(c_stable) 0.6];
    plot([ci_L(i) ci_H(i)], [yy yy], '-', 'Color', line_col, 'LineWidth', lw);
    plot([ci_L(i) ci_L(i)], [yy-cap yy+cap], '-', 'Color', line_col, 'LineWidth', lw);
    plot([ci_H(i) ci_H(i)], [yy-cap yy+cap], '-', 'Color', line_col, 'LineWidth', lw);
end
scatter(x_stable, y_list, 20, 'filled', 'o', 'MarkerFaceColor', c_stable, 'MarkerEdgeColor', 'none');
scatter(x_move, y_list, 20, 'filled', 'o', 'MarkerFaceColor', c_move, 'MarkerEdgeColor', 'none');
end

function rgb = hex2rgb(hexStr)
hexStr = regexprep(hexStr, '#', '');
rgb = sscanf(hexStr, '%2x%2x%2x', [1 3]) / 255;
end
