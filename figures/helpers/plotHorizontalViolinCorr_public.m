function plotHorizontalViolinCorr_public(ax, null_dist, real_val, p_val, y_pos, color_violin)

if nargin < 6
    color_violin = [217 227 243]/255;
end

c_Sig = [171, 77, 62] / 255; 
c_NotSig = [113, 129, 158] / 255;

if p_val < 0.05
    current_color = c_Sig;
else
    current_color = c_NotSig;
end

violin_width = 0.35;
marker_size = 10;

hold(ax, 'on');

[f, xi] = ksdensity(null_dist);
f_normalized = f / max(f) * violin_width;

x_coords = [xi, fliplr(xi)];
y_coords = [y_pos + f_normalized, y_pos - fliplr(f_normalized)];

patch(ax, x_coords, y_coords, color_violin, ...
      'FaceAlpha', 1, 'EdgeColor', 'none');

median_null = median(null_dist);
plot(ax, [median_null, median_null], [y_pos - violin_width, y_pos + violin_width], ...
     'Color', [0.2 0.2 0.2], 'LineWidth', 0.85);

plot(ax, real_val, y_pos, 'o', ...
    'MarkerFaceColor', current_color, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', marker_size);


hold(ax, 'off');

end
