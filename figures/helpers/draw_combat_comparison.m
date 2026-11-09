function fig = draw_combat_comparison(x, y, show_ref_line)
    figure('Name', 'ComBat_Comparison', 'Position', [100  100  340  340], 'Color', 'white');    
    
    if nargin < 3
        show_ref_line = false;
    end
    
    valid_idx = ~isnan(x) & ~isnan(y);
    x_clean = x(valid_idx);
    y_clean = y(valid_idx);
    
    global_max = max(abs([x_clean; y_clean]));
    limit_val = global_max * 1.1; 
    axis_lims = [-limit_val, limit_val];
    
    if show_ref_line
        plot(axis_lims, axis_lims, 'r--', 'LineWidth', 0.85);
        hold on;
    end
    
    scatter(x, y, 15, [128,115,172]/255, 'filled', 'MarkerFaceAlpha', 0.8);
    hold on;
    
    p_poly = polyfit(x_clean, y_clean, 1);
    x_fit = linspace(axis_lims(1), axis_lims(2), 100);
    y_fit = polyval(p_poly, x_fit);
    plot(x_fit, y_fit, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.85);
    
    xlim(axis_lims);
    ylim(axis_lims);
    
    [r, p_val] = corr(x_clean, y_clean, 'Rows', 'complete');
    
    txt_color = 'k';
    if abs(r) > 0.9, txt_color = [0 0.5 0]; end
    if abs(r) < 0.1, txt_color = [0.5 0.5 0.5]; end
    
    if p_val < 0.001
        p_str = 'p < 0.001';
    else
        p_str = sprintf('p = %.3e', p_val);
    end
    
    label_str = {sprintf('r = %.3f', r), p_str};
    text(axis_lims(1) + 0.05*range(axis_lims), ...
        axis_lims(2) - 0.10*range(axis_lims), ...
        label_str, ...
        'Color', txt_color, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');
        
    grid off; box off;
    set(gca, 'TickDir', 'out', 'LineWidth', 0.85, 'FontSize', 8, 'TickLength', 4*[0.01 0.025]);
    axis square;
    fig = gcf;
end

