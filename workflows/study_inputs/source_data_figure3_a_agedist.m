function source_data_figure3_a_agedist(brain_file,combat_dir,history_file,out)
if ~isfolder(out),mkdir(out);end
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','figures','helpers'));

string_dem = {'area','thickness','volume'};
dem = 2;
load(brain_file,'brain_data');
brain_data.pc7_2 = -brain_data.pc7_2;
brain_data.pc7_2_recruit = -brain_data.pc7_2_recruit;

all_vars = fieldnames(brain_data);
area_var = all_vars(contains(all_vars,'area') & ~contains(all_vars,'500_aparc_'));
areadata = nan(length(brain_data.age),length(area_var));
for arvv = 1:length(area_var)
    areadata(:,arvv) = brain_data.(area_var{arvv});
end
brain_data.Total_area = sum(areadata,2);

brain_data=attach_move_time(brain_data,history_file);

savepathss = [combat_dir filesep];
Culstr = '_cul'; SESstr = ''; Envstr = '_withoutEnv';
load(fullfile(combat_dir,sprintf(['data_harmonized_' string_dem{dem} '_demo%s_valididx%s%s.mat'], Culstr, SESstr, Envstr)), 'valididx_base');
idx_use = valididx_base;

LivedPlace = brain_data.livedplace_32(idx_use);
Age_Current = brain_data.age(idx_use);

is_stayer = (LivedPlace == 1);
is_mover  = (LivedPlace >= 2);

TimeSince = brain_data.TimeSinceMove(idx_use);
Age_Move_All = Age_Current(is_mover) - TimeSince(is_mover);
valid_mask = Age_Move_All > 0 & Age_Move_All <= Age_Current(is_mover);
MoveAge_Data = Age_Move_All(valid_mask);

Fig3a = struct();

[counts_move, edges_move] = histcounts(MoveAge_Data, 'BinMethod', 'integers', 'Normalization', 'probability');
centers_move = edges_move(1:end-1) + 0.5;
cdf_move = cumsum(counts_move);

Fig3a.Plot1_centers = centers_move;
Fig3a.Plot1_counts = counts_move;
Fig3a.Plot1_cdf = cdf_move;

bin_width = 2;
edges_age = floor(min(Age_Current)):bin_width:ceil(max(Age_Current));
centers_age = (edges_age(1:end-1) + edges_age(2:end)) / 2;

counts_stayer = histcounts(Age_Current(is_stayer), edges_age, 'Normalization', 'probability');
counts_mover  = histcounts(Age_Current(is_mover), edges_age, 'Normalization', 'probability');

[f_s, xi_s] = ksdensity(Age_Current(is_stayer));
[f_m, xi_m] = ksdensity(Age_Current(is_mover));

N_stayer = sum(is_stayer);
N_mover  = sum(is_mover);
mean_age_stayer = mean(Age_Current(is_stayer));
mean_age_mover  = mean(Age_Current(is_mover));
[~, p_ranksum] = ranksum(Age_Current(is_stayer), Age_Current(is_mover));

Fig3a.Plot2_centers = centers_age;
Fig3a.Plot2_counts_stayer = counts_stayer;
Fig3a.Plot2_counts_mover = counts_mover;
Fig3a.Plot2_bin_width = bin_width;
Fig3a.Plot2_kde_xi_stayer = xi_s;
Fig3a.Plot2_kde_f_stayer = f_s;
Fig3a.Plot2_kde_xi_mover = xi_m;
Fig3a.Plot2_kde_f_mover = f_m;
Fig3a.Plot2_N_stayer = N_stayer;
Fig3a.Plot2_N_mover = N_mover;
Fig3a.Plot2_mean_age_stayer = mean_age_stayer;
Fig3a.Plot2_mean_age_mover = mean_age_mover;
Fig3a.Plot2_p_ranksum = p_ranksum;

colorpair = {'#f4a582','#E6A094','#d6604d';   % Red
             '#92c5de','#8DA1C6','#4163a0'};   % Blue
Fig3a.colorpair = colorpair;

save_path = fullfile(out,'Fig3a_AgeDist_SourceData.mat');
save(save_path, '-struct', 'Fig3a');
fprintf('Source Data for Figure 3a saved to:\n  %s\n', save_path);

end

