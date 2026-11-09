function splot_hemis_bound(X, Y, ColorScatter, ColorEdge, ColorLine, Colorbound, lhidx, p_spin, n_rot, target_xlim, custom_stats_text)

    math_minus = char(8722);

    scatter(X(lhidx), Y(lhidx), 'MarkerFaceColor', ColorScatter{1}, 'MarkerEdgeColor', ColorEdge, ...
        'MarkerFaceAlpha', 1, 'SizeData', 15);
    hold on;
    scatter(X(~lhidx), Y(~lhidx), 'MarkerFaceColor', ColorScatter{2}, 'MarkerEdgeColor', ColorEdge, ...
        'MarkerFaceAlpha', 1, 'SizeData', 15);
    
    X = X(:); Y = Y(:);
    
    [p_fit, S_fit] = polyfit(X, Y, 1);
    
    if nargin > 9 && ~isempty(target_xlim)
        xlim(target_xlim);
    end

    x_range = xlim;

    x_fit = linspace(x_range(1), x_range(2), 100)';
    
    [y_fit, delta] = polyconf(p_fit, x_fit, S_fit, 'PredOpt', 'curve', 'SimOpt', 'on');
    
    plot(x_fit, y_fit, 'Color', ColorLine, 'LineWidth', 2);
    fill([x_fit; flip(x_fit)], [y_fit - delta; flip(y_fit + delta)], Colorbound, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    valid_idx = ~isnan(X) & ~isnan(Y);
    [r, p_pearson] = corr(X(valid_idx), Y(valid_idx));
    
    
    use_spin_p = false;
    if nargin > 7 && ~isempty(p_spin)
        p_show = p_spin;
        use_spin_p = true;
    else
        p_show = p_pearson;
    end
    
    if nargin < 9 || isempty(n_rot)
        n_rot = 1000;
    end
    
    r_str = sprintf('%.3f', r);
    r_str = strrep(r_str, '-', math_minus);
    
    if use_spin_p && p_show == 0
        p_str = sprintf('< %.3f', 1/n_rot);
        p_label = 'p_{spin}';
    else
        if p_show >= 0.0001
            p_str = sprintf('= %.3f', p_show);
        else
            exponent = floor(log10(p_show));
            mantissa = p_show / 10^exponent;
            p_str_sci = sprintf('%.1f \\times 10^{%d}', mantissa, exponent);
            p_str = ['= ' strrep(p_str_sci, '-', math_minus)];
        end
        
        if use_spin_p
            p_label = 'p_{spin}';
        else
            p_label = 'p';
        end
    end
    
    if nargin > 10 && ~isempty(custom_stats_text)
        stats_text = custom_stats_text;
    else
        stats_text = sprintf('\\it{r} = %s, \\it{%s} %s', r_str, p_label, p_str);
    end
    
    text(0.90, 0.92, stats_text, ...
         'Units', 'normalized', ...
         'FontSize', 7, ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', ...
         'Interpreter', 'tex', ...
         'BackgroundColor', 'none');

    ax = gca;
    ax.FontSize = 8;
    ax.FontName = 'Arial';
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.LineWidth = 0.85;

    hold off;
    grid off; axis square;
    set(gca, 'TickDir', 'out', 'LineWidth', 0.85, 'FontSize', 8);
    set(gca, 'TickLength', 4*[0.01 0.025]);

end
