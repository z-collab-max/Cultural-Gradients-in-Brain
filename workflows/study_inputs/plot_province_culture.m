function plot_province_culture(source,geometry,out)
if isempty(geometry),geometry=fullfile(fileparts(mfilename('fullpath')),'..','..','assets','geography','china_shapes');end
if ~isfolder(out),mkdir(out);end
set(groot,'defaultFigureVisible','off','defaultAxesFontName','Arial','defaultAxesFontSize',8,'defaultAxesLineWidth',.85);
load(source);
shp_path_poly=fullfile(geometry,'bou2_4p_utf8.shp');
shp_path_line=fullfile(geometry,'bou1_4l.shp');
titles_map = {'Origin Culture (Stayers)', 'Destination Culture (Relocators)'};
names_map = {'Stayers', 'Relocators'};
for ccc = 1:2
    Prov_Status   = eval(sprintf('PC%d_Prov_Status', ccc));
    global_median = eval(sprintf('PC%d_global_median', ccc));
    clim_range_pc = eval(sprintf('PC%d_clim_range', ccc));
    for map_idx = 1:2
        figure('Name', sprintf('ChinaMap_Culture_PC%d_%s', ccc, names_map{map_idx}), 'Color', 'w', 'Position', [100 100 800 600]);
        ax = worldmap([18 54], [73 136]);
        setm(ax, 'MapProjection', 'eqaconic');
        
        provinces = shaperead(shp_path_poly, 'UseGeoCoords', true);
        bordersl  = shaperead(shp_path_line, 'UseGeoCoords', true);
        geoshow(bordersl, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.85);
        
        face_cols = nan(numel(provinces), 3);
        
        for k = 1:length(provinces)
            raw_name = provinces(k).NAME;
            pid_found = 0;
            for pid = 1:34
                short_name = Prov_List{pid}(1:2);
                if contains(raw_name, short_name)
                    pid_found = pid;
                    break;
                end
            end
            this_color = [0.95 0.95 0.95]; % default light gray
            if pid_found > 0
                val = Prov_Status(pid_found, map_idx);
                
                if ~isnan(val)
                    if val >= global_median
                        raw_ratio = (val - global_median) / (clim_range_pc(2) - global_median);
                        raw_ratio = min(max(raw_ratio, 0), 1);
                        adj_ratio = start_sat + (1 - start_sat) * raw_ratio;
                        this_color = (1 - adj_ratio) * color_Mid + adj_ratio * color_High_End;
                    else
                        raw_ratio = (global_median - val) / (global_median - clim_range_pc(1));
                        raw_ratio = min(max(raw_ratio, 0), 1);
                        adj_ratio = start_sat + (1 - start_sat) * raw_ratio;
                        this_color = (1 - adj_ratio) * color_Mid + adj_ratio * color_Low_End;
                    end
                end
            end
            face_cols(k, :) = this_color;
        end
        
        faceColors = makesymbolspec('Polygon', ...
            {'INDEX', [1 numel(provinces)], ...
            'FaceColor', face_cols, ...
            'EdgeColor', [1 1 1], ...
            'LineWidth', 0.85});
        geoshow(provinces, 'DisplayType', 'polygon', 'SymbolSpec', faceColors);
        
        setm(gca, 'Frame', 'off', 'Grid', 'off', 'ParallelLabel', 'off', ...
            'MeridianLabel', 'off', 'FFaceColor', 'w');
        axis off;
        title(sprintf('Culture PC%d — %s', ccc, titles_map{map_idx}), 'FontSize', 8, 'FontName', 'Arial');
        n_steps = 100;
        ratio_vec = linspace(start_sat, 1, n_steps)';
        blues = (1 - ratio_vec) .* color_Mid + ratio_vec .* color_Low_End;
        blues = flipud(blues);
        reds = (1 - ratio_vec) .* color_Mid + ratio_vec .* color_High_End;
        
        custom_cmap = [blues; reds];
        
        c = colorbar(ax, 'Position', [0.85 0.3 0.02 0.4]);
        colormap(ax, custom_cmap);
        caxis(ax, clim_range_pc);
        c.Label.String = 'Culture Score';
        set(c, 'Box', 'off');
    end
end

figs=findall(groot,'Type','figure');
for j=1:numel(figs)
 set(figs(j),'Position',[100 100 1050 650]);
 name=regexprep(get(figs(j),'Name'),'[^a-zA-Z0-9_]','_');if isempty(name),name=sprintf('panel%d',j);end
 exportgraphics(figs(j),fullfile(out,['plot_fig3h_' name '.png']),'Resolution',150);
 savefig(figs(j),fullfile(out,['plot_fig3h_' name '.fig']));
end
close all;
end
