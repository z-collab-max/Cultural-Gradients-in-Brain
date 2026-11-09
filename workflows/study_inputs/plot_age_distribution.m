function plot_age_distribution(source,out)
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'..','..','figures','helpers'),fullfile(here,'..','..','simulation'));oldfig=findall(groot,'Type','figure');
load(source);

% --- Color setup ---
rgb_1 = validatecolor(colorpair{1,3}); Color_Line1 = rgb_1 * 0.8;
rgb_2 = validatecolor(colorpair{2,3}); Color_Line2 = rgb_2 * 0.8;
Color_Stayer = Color_Line1;
Color_Mover  = Color_Line2;

%% ========================================================================
% Plot 1: Distribution of Age at First Relocation
% =========================================================================
figure('Name', 'Age_at_First_Relocation', 'Color', 'w', 'Position', [100, 250, 0.9*0.8*0.8*650, 150]);

yyaxis left
% Proportion area + line
h_area = area(Plot1_centers, Plot1_counts, 'FaceColor', Color_Line1, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
hold on;
plot(Plot1_centers, Plot1_counts, 'Color', Color_Line1, 'LineWidth', 0.85);
ylabel('Proportion', 'FontSize', 8, 'FontWeight', 'bold');
set(gca, 'YColor', Color_Line1);
ylim([0 max(Plot1_counts)*1.2]);

yyaxis right
% Cumulative line
plot(Plot1_centers, Plot1_cdf, '-', 'Color', Color_Line2, 'LineWidth', 0.85);
ylabel('Cumulative proportion', 'FontSize', 8);
set(gca, 'YColor', Color_Line2);
ylim([0 1.05]);

% Annotations
title('Distribution of age at first relocation', 'FontSize', 11);
xlabel('Age (years)', 'FontSize', 8);
box off;

% Mark key ages (19 and 25)
for target_age = [19, 25]
    [~, idx] = min(abs(Plot1_centers - target_age));
    if abs(Plot1_centers(idx) - target_age) < 0.5
        val_cdf = Plot1_cdf(idx);
        text(target_age, val_cdf, sprintf(' Age %d\n CDF: %.1f%%', target_age, val_cdf*100), ...
            'Color', Color_Line2, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontSize', 7);
        xline(target_age, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.85, 'HandleVisibility', 'off');
    end
end

set(gca, 'XLim', [0 35], 'XTick', 0:5:35);
box off;
try nice_axis(gca); catch; end

%% ========================================================================
% Plot 2: Current Age Distribution — Stayers vs. Movers (Grouped Bar)
% =========================================================================
figure('Name', 'Current_Age_Stayers_vs_Movers', 'Color', 'w', 'Position', [100, 250, 0.9*0.8*0.8*650, 150]);
hold on;

% Grouped bar chart
bar_data = [Plot2_counts_stayer(:), Plot2_counts_mover(:)];
b = bar(Plot2_centers, bar_data, 'grouped');

% Stayers style
b(1).FaceColor = Color_Stayer;
b(1).FaceAlpha = 0.3;
b(1).EdgeColor = 'w';
b(1).LineWidth = 0.5;
b(1).DisplayName = sprintf('Stayers (N=%d)', Plot2_N_stayer);

% Movers style
b(2).FaceColor = Color_Mover;
b(2).FaceAlpha = 0.3;
b(2).EdgeColor = 'w';
b(2).LineWidth = 0.5;
b(2).DisplayName = sprintf('Movers (N=%d)', Plot2_N_mover);

% Kernel density curves
plot(Plot2_kde_xi_stayer, Plot2_kde_f_stayer * Plot2_bin_width, ...
    'Color', Color_Stayer, 'LineWidth', 0.85, 'HandleVisibility', 'off');
plot(Plot2_kde_xi_mover, Plot2_kde_f_mover * Plot2_bin_width, ...
    'Color', Color_Mover, 'LineWidth', 0.85, 'HandleVisibility', 'off');

% Decorations
title('Current Age Distribution: Stayers vs. Movers', 'FontSize', 14);
xlabel('Current Age (years)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Proportion', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'Box', 'off');
xlim([0 32]);
try nice_axis(gca); catch; end

if ~isfolder(out),mkdir(out);end
figs=setdiff(findall(groot,'Type','figure'),oldfig);
for j=1:numel(figs)
 set(figs(j),'Position',[100 100 1050 650]);
 name=regexprep(get(figs(j),'Name'),'[^a-zA-Z0-9_]','_');if isempty(name),name=sprintf('panel%d',j);end
 publication_style(figs(j));exportgraphics(figs(j),fullfile(out,['plot_fig3a_' name '.png']),'Resolution',150);
 savefig(figs(j),fullfile(out,['plot_fig3a_' name '.fig']));
end
close(figs);
end


