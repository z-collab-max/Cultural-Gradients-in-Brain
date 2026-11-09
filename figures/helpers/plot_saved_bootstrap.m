function [h, p_val] = plot_saved_bootstrap(null_dist, real_val, title_str, xlabel_str, saved_p, digits)
if nargin<6,digits=4;end
if nargin < 5
    saved_p = NaN;
end
hold on;

h = histogram(null_dist, 30, 'Normalization', 'pdf', ...
    'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Null Dist.');

[f_kde, xi_kde] = ksdensity(null_dist);
plot(xi_kde, f_kde, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.85, 'DisplayName', 'KDE Fit');

xlim([-1, 1]);
xticks([-1, -0.5, 0, 0.5, 1]);

current_ylim = get(gca, 'YLim');
max_y = current_ylim(2) * 1.05;
ylim([0, max_y]);

p_val = min(1,2*min(mean(null_dist<=real_val),mean(null_dist>=real_val)));
if p_val < 0.05
    line_style = '-';
    line_width = 1.75;
else
    line_style = '--';
    line_width = 1.5;
end
plot([real_val real_val], [0, max_y], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', line_width, 'LineStyle', line_style);

if p_val == 0
    p_str = sprintf('\\it{p_{boot}} < %.4f', 1/length(null_dist));
else
    p_str = sprintf('\\it{p_{boot}} = %.*f', digits, p_val);
end


if iscell(title_str)
    final_title = title_str;
else
    final_title = {title_str;};
end
stats_text = sprintf('%s', p_str);


text(0.6, 0.92, stats_text, ...
    'Units', 'normalized', ...
    'FontSize', 7, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Interpreter', 'tex', ...
    'BackgroundColor', 'none');


title(['\fontsize{10.5}\rm' final_title]);
xlabel(xlabel_str, 'FontSize', 8);
ylabel('Probability Density', 'FontSize', 8);




box off;
axis square;
set(gca,'TickDir','out','LineWidth',.85,'FontSize',8,'FontName','Arial');
end
