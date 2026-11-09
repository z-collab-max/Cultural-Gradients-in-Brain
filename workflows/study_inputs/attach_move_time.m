function brain_data=attach_move_time(brain_data,filename)


if ~isfile(filename)
    error('文件 %s 未找到，请确认路径是否正确。', filename);
end

fprintf('正在读取数据文件: %s ...\n', filename);

opts = detectImportOptions(filename);

raw_data = readtable(filename, opts);

required_vars = {'id', 'year', 'PROVINCE', 'CITY'}; 
raw_data.Properties.VariableNames = upper(raw_data.Properties.VariableNames); 

if ~all(ismember(upper(required_vars), raw_data.Properties.VariableNames))
    error('CSV文件缺少必要的列 (id, year, PROVINCE, CITY)');
end

unique_ids = unique(raw_data.ID);
n_subs = length(unique_ids);

fprintf('共检测到 %d 个被试，开始计算...\n', n_subs);

out_ID = unique_ids;
out_Age = nan(n_subs, 1);
out_IsMover = false(n_subs, 1);
out_TimeSinceMove = nan(n_subs, 1); % 核心变量
out_BirthProv = cell(n_subs, 1);
out_ScanProv = cell(n_subs, 1);

wait_h = waitbar(0, '正在处理居住史...');

for i = 1:n_subs
    curr_id = unique_ids(i);
    
    idx_rows = raw_data.ID == curr_id;
    sub_table = raw_data(idx_rows, :);
    
    [~, sort_idx] = sort(sub_table.YEAR);
    sub_table = sub_table(sort_idx, :);
    
    years = sub_table.YEAR;
    provinces = sub_table.PROVINCE; % 使用省份作为文化单位 (如需城市请改用 CITY)
    
    birth_year = years(1);
    scan_year = years(end);
    current_age = scan_year - birth_year; % 近似扫描年龄
    
    birth_place = char(provinces{1});   % 出生地 (转为char方便比较)
    current_place = char(provinces{end}); % 扫描地
    
    out_Age(i) = current_age;
    out_BirthProv{i} = birth_place;
    out_ScanProv{i} = current_place;
    
    
    
    diff_idx = find(~strcmp(provinces, birth_place), 1, 'first');
    
    if ~isempty(diff_idx)
        out_IsMover(i) = true;
        
        year_first_move = years(diff_idx);
        
        duration = scan_year - year_first_move;
        
        if duration < 0, duration = 0; end
        
        out_TimeSinceMove(i) = duration;
    else
        out_IsMover(i) = false;
        out_TimeSinceMove(i) = 0;
    end
    
    if mod(i, 500) == 0
        waitbar(i/n_subs, wait_h, sprintf('已处理 %d / %d', i, n_subs));
    end
end
close(wait_h);

ResultTable = table(out_ID, out_Age, out_IsMover, out_TimeSinceMove, out_BirthProv, out_ScanProv, ...
    'VariableNames', {'id', 'Age_Calc', 'Is_Mover', 'Time_Since_Move', 'Birth_Province', 'Scan_Province'});

fprintf('\n=== 计算完成 ===\n');
fprintf('Movers (离开过出生地): %d\n', sum(out_IsMover));
fprintf('Stayers (始终在出生地): %d\n', sum(~out_IsMover));
fprintf('Movers 平均离开时长: %.2f 年\n', mean(out_TimeSinceMove(out_IsMover)));



brain_data.TimeSinceMove = nan(length(brain_data.MRIID), 1);
brain_data.IsMover_Calc = nan(length(brain_data.MRIID), 1);

[Lia, Locb] = ismember(brain_data.MRIID, ResultTable.id);

brain_data.TimeSinceMove(Lia) = ResultTable.Time_Since_Move(Locb(Lia));
brain_data.IsMover_Calc(Lia) = double(ResultTable.Is_Mover(Locb(Lia)));

fprintf('已将计算结果合并入 brain_data.\n');
end
