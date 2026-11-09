
function plot_ed2(input_file)
out_dir = fileparts(mfilename('fullpath'));

package=fileparts(fileparts(mfilename('fullpath')));addpath(fileparts(mfilename('fullpath')));
out_dir=fullfile(package,'outputs','simulation_figures','ed2');if ~exist(out_dir,'dir'),mkdir(out_dir);end
suffix='';
if nargin<1,input_file=fullfile(package,'simulation','reference','ed2','EDFig2c_simulation.mat');end
S=load(input_file);
nBins = numel(S.artif_rhos);

colors = [107 152 187;
          105 170 162;
          151 143 184;
          209 168 102;
          217 135 124] / 255;
ref_color = [184 189 193] / 255;
map_size = 7;
map_alpha = 0.58;
cal_size = 7;
cal_alpha = 0.62;
fit_line_width = 0.85;
ref_line_width = 0.60;
fit_length_factor = 1.25;
map_lim = [-4, 4];
map_ticks = -4:2:4;
cal_lim = [-1.1, 1.1];
cal_ticks = [-1 -0.5 0 0.5 1];

fig = figure('Color', 'w', 'Position', [40 250 2700 560]);
panel_x = [0.025 0.120 0.225 0.320 0.465 0.560 0.665 0.760];
panel_w = 0.085;
panel_y = 0.19;
panel_h = 0.66;

arms = {'M', 'C'};
short = {'s', 'e'};
stages = {'start', 'end'};
stage_titles = {'Raw', 'De-artifacted'};
example_isim_all = cell(2, 1);
legend_handles = gobjects(nBins, 1);

for a = 1:2
    r_raw = S.(['corr_est_start_GT_' arms{a}]);
    med_raw = median(r_raw, 2);
    [~, example_isim] = min(abs(r_raw - med_raw), [], 2);
    example_isim_all{a} = example_isim;

    for local_panel = 1:4
        global_panel = (a - 1) * 4 + local_panel;
        ax = axes('Parent', fig, 'Position', ...
            [panel_x(global_panel), panel_y, panel_w, panel_h]);
        hold(ax, 'on');

        if local_panel <= 2
            st = local_panel;
            plot(ax, map_lim, map_lim, '--', 'Color', ref_color, ...
                'LineWidth', ref_line_width);

            map_est = S.(['maps_' short{st} arms{a}]);
            fit_coef = nan(nBins, 2);
            x_show = cell(nBins, 1);
            y_show = cell(nBins, 1);
            for ibin = 1:nBins
                row_idx = (example_isim(ibin) - 1) * nBins + ibin;
                x_show{ibin} = zscore(S.maps_GT(row_idx, :)');
                y_show{ibin} = zscore(map_est(row_idx, :)');
                fit_coef(ibin, :) = polyfit(x_show{ibin}, y_show{ibin}, 1);
            end

            fit_x = [min(cellfun(@min, x_show)), max(cellfun(@max, x_show))];
            fit_mid = mean(fit_x);
            fit_half = diff(fit_x) * fit_length_factor / 2;
            fit_x = fit_mid + [-fit_half, fit_half];
            for ibin = 1:nBins
                plot(ax, fit_x, polyval(fit_coef(ibin, :), fit_x), '-', ...
                    'Color', colors(ibin, :), 'LineWidth', fit_line_width);
            end
            for ibin = 1:nBins
                scatter(ax, x_show{ibin}, y_show{ibin}, map_size, colors(ibin, :), ...
                    'filled', 'MarkerFaceAlpha', map_alpha, 'MarkerEdgeColor', 'none');
            end

            xlim(ax, map_lim); ylim(ax, map_lim);
            set(ax, 'XTick', map_ticks, 'YTick', map_ticks);
            xlabel(ax, 'ground truth coefficient (z-score)');
            if local_panel == 1
                ylabel(ax, {'estimated coefficient', '(z-score)'});
            end
            title(ax, stage_titles{st});
        else
            st = local_panel - 2;
            data_y = S.(['corr_est_' stages{st} '_SA_' arms{a}]);
            plot(ax, cal_lim, cal_lim, '--', 'Color', ref_color, ...
                'LineWidth', ref_line_width);
            plot(ax, cal_lim, [0 0], ':', 'Color', ref_color, ...
                'LineWidth', ref_line_width);
            plot(ax, [0 0], cal_lim, ':', 'Color', ref_color, ...
                'LineWidth', ref_line_width);

            fit_x = [min(S.corr_GT_SA(:)), max(S.corr_GT_SA(:))];
            fit_mid = mean(fit_x);
            fit_half = diff(fit_x) * fit_length_factor / 2;
            fit_x = fit_mid + [-fit_half, fit_half];
            for ibin = 1:nBins
                fit_coef = polyfit(S.corr_GT_SA(ibin, :), data_y(ibin, :), 1);
                plot(ax, fit_x, polyval(fit_coef, fit_x), '-', ...
                    'Color', colors(ibin, :), 'LineWidth', fit_line_width);
            end
            for ibin = 1:nBins
                h = scatter(ax, S.corr_GT_SA(ibin, :), data_y(ibin, :), ...
                    cal_size, colors(ibin, :), 'filled', ...
                    'MarkerFaceAlpha', cal_alpha, 'MarkerEdgeColor', 'none');
                if a == 2 && st == 2
                    legend_handles(ibin) = h;
                end
            end

            xlim(ax, cal_lim); ylim(ax, cal_lim);
            set(ax, 'XTick', cal_ticks, 'YTick', cal_ticks);
            xlabel(ax, 'r (cultural gradient, SA)');
            if local_panel == 3
                ylabel(ax, 'r (cultural gradient, SA)');
            end
            title(ax, stage_titles{st});
        end

        axis(ax, 'square');
        set(ax, 'TickDir', 'out', 'LineWidth', 0.5, 'FontSize', 7, ...
            'TickLength', 2 * [0.01 0.025]);
        xtickangle(ax, 0);
    end
end


r_sM = S.corr_est_start_GT_M; r_eM = S.corr_est_end_GT_M;
r_sC = S.corr_est_start_GT_C; r_eC = S.corr_est_end_GT_C;
exM = example_isim_all{1}; exC = example_isim_all{2};
legend_text = cell(nBins, 1);
for ibin = 1:nBins
    legend_text{ibin} = sprintf('%.1f   %.2f/%.2f   %.2f/%.2f', ...
        S.artif_rhos(ibin), ...
        r_sM(ibin, exM(ibin)), r_eM(ibin, exM(ibin)), ...
        r_sC(ibin, exC(ibin)), r_eC(ibin, exC(ibin)));
end

out_base = fullfile(out_dir, 'EDFig2c_simulation');
publication_style(fig);savefig(fig, [out_base '.fig']);
saveas(fig, [out_base '.svg']);
exportgraphics(fig, [out_base '.png'], 'Resolution', 300);
fprintf('Saved %s.svg and %s.png\n', out_base, out_base);
close(fig);
end

