function source_data_figure3_e_chinamap(brain_file,combat_dir,history_file,out)
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

controlvars = {'sex','age','edu','Gpc1','Gpc2','Gpc3','Gpc4','Gpc5','Gpc6','Gpc7','Gpc8','Gpc9','Gpc10'};

savepathss = [combat_dir filesep];
Culstr = '_cul'; SESstr = ''; Envstr = '_withoutEnv';
load(fullfile(combat_dir,sprintf(['data_harmonized_' string_dem{dem} '_demo%s_valididx%s%s.mat'], Culstr, SESstr, Envstr)), 'valididx_base');
idx_use = valididx_base;

Usenames_pc = {'pc7_1','pc7_2'};

Prov_List = {'北京市','天津市','河北省','山西省','内蒙古自治区',...
             '辽宁省','吉林省','黑龙江省','上海市','江苏省',...
             '浙江省','安徽省','福建省','江西省','山东省',...
             '河南省','湖北省','湖南省','广东省','广西壮族自治区',...
             '海南省','重庆市','四川省','贵州省','云南省',...
             '西藏自治区','陕西省','甘肃省','青海省','宁夏回族自治区',...
             '新疆维吾尔自治区','台湾省','香港特别行政区','澳门特别行政区'};

idx_stay = idx_use & (brain_data.province_recruit == brain_data.province_birth);
idx_move = idx_use & (brain_data.livedplace_32 == 2);

Fig3e = struct();
Fig3e.Prov_List = Prov_List;

for ccc = 1:2
    n_provs = 34;
    Prov_Status = nan(n_provs, 2); % Col1=Stayer Origin, Col2=Mover Destination
    
    for pid = 1:34
        idx_born_in_p = (brain_data.province_birth == Prov_List{pid});
        
        data_stay = brain_data.(Usenames_pc{ccc})(idx_stay & idx_born_in_p);
        if ~isempty(data_stay)
            Prov_Status(pid, 1) = mean(data_stay, 'omitnan');
        end
        
        data_move = brain_data.([Usenames_pc{ccc} '_recruit'])(idx_move & idx_born_in_p);
        if ~isempty(data_move)
            Prov_Status(pid, 2) = mean(data_move, 'omitnan');
        end
    end
    
    assert(any(isfinite(Prov_Status(:,1))) && any(isfinite(Prov_Status(:,2))), 'No provinces matched');
    all_vals = Prov_Status(~isnan(Prov_Status));
    global_median = median(brain_data.(Usenames_pc{ccc})(idx_use));
    max_dev = max(abs(all_vals - global_median));
    clim_range = [global_median - max_dev, global_median + max_dev];
    
    Fig3e.(sprintf('PC%d_Prov_Status', ccc)) = Prov_Status;
    Fig3e.(sprintf('PC%d_global_median', ccc)) = global_median;
    Fig3e.(sprintf('PC%d_clim_range', ccc)) = clim_range;
end

Fig3e.color_High_End = hex2rgb('#E6A094');
Fig3e.color_Low_End  = hex2rgb('#8DA1C6');
Fig3e.color_Mid      = [1 1 1];
Fig3e.start_sat      = 0.25;

save_path = fullfile(out,'Fig3e_ChinaMap_SourceData.mat');
save(save_path, '-struct', 'Fig3e');
fprintf('Source Data for Figure 3e saved to:\n  %s\n', save_path);

end

