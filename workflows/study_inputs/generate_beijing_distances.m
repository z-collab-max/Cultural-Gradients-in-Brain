function generate_beijing_distances(culture_file,extra_file,out,draws,seed)
if nargin<4,draws=10000;end
if nargin<5,seed=42;end
if ~isfolder(out),mkdir(out);end
rng(seed,'twister');
fzscore = @(x) (x - mean(x))/std(x);
load(culture_file);

%%
tp = readtable(extra_file,'Sheet','Sheet1');
tpdata1 = table2struct(tp,'ToScalar',true);
tp = readtable(extra_file,'Sheet','Sheet2');
tpdata2 = table2struct(tp,'ToScalar',true);

tpdata = struct;
vnames = fieldnames(tpdata1);
for vvv = 1:length(vnames)
    tpdata.(vnames{vvv}) = [tpdata1.(vnames{vvv});tpdata2.(vnames{vvv})];
end

XJidx = tpdata.province == 31;
tpdata.provincei(XJidx) = {'新疆维吾尔自治区'};
TBidx = tpdata.province == 26;
tpdata.provincei(TBidx) = {'西藏自治区'};
tpdata.province = tpdata.provincei;
tpdata.province = categorical(tpdata.province);
tpdata.province = tpdata.province(:);

tpdata.ocity_y = str2double(tpdata.ocity_y);
tpdata.age = str2double(tpdata.age);
%%
vnames = fieldnames(tpdata);
Culdata_cat = struct;
for vvv = 1:length(vnames)-1
    Culdata_cat.(vnames{vvv}) = [Culdata.(vnames{vvv});tpdata.(vnames{vvv})];
end

%% 剔除非汉族被试
Culdata_cal_han = struct;
data = struct;
tmpfields = fieldnames(Culdata_cat);
hanidx = Culdata_cat.race_han_yn == 1;
localidx = ((Culdata_cat.city_fa == 1)&(Culdata_cat.city_ma == 1))|Culdata_cat.province == '新疆维吾尔自治区'|Culdata_cat.province == '西藏自治区';
%     transidx = (isnan(Culdata_cat.ocity_where))|((~isnan(Culdata_cat.ocity_where)&(Culdata_cat.ocity_y<=5)));

for varid = 1:length(tmpfields)
    data.(tmpfields{varid}) = Culdata_cat.(tmpfields{varid})(hanidx&localidx);
end

%% 计算其他的控制变量
%%%%%%%%%%%%%% 筛选变量 %%%%%%%%%%%%%
%%% 民族（1:汉）
Culdata_cal_han.race_han_yn = (data.race_han_yn);
%%% 宗教（2：无）
Culdata_cal_han.religion_yn = (data.religion_yn);
%%% 上大学前生活的城市类别（1：直辖市或者省会城市 2：中小城市 3：县城 4：乡镇一级城市）
Culdata_cal_han.city_beforecoll = (data.city_beforecoll);
%%% 在其他省生活过吗 生活过多少年
Culdata_cal_han.ocity_yn = (data.ocity_yn);
Culdata_cal_han.ocity_where = (data.ocity_where);
Culdata_cal_han.ocity_y = (data.ocity_y);

%%%%%%%%%%%%%% 控制变量 %%%%%%%%%%%%%%


%%% 省份
Culdata_cal_han.province = data.province;
%%% 城市
%     Culdata_cal_han.city = data.city;
%%% 年龄
Culdata_cal_han.age = (data.age);
%%% 性别
Culdata_cal_han.gender = (data.gender);
%%% 专业
Culdata_cal_han.majorincoll = (data.majorincoll);
%%% 受教育年限 1 <=10 2 11 3 12 4 13 5 14 6 15 7 16 8 17 9 >=18
Culdata_cal_han.edu_year = data.edu_year;

%%% 父母受教育程度
Culdata_cal_han.edu_ma = (data.edu_ma);
Culdata_cal_han.edu_fa = (data.edu_fa);
%%% 父母城市
Culdata_cal_han.city_ma = (data.city_ma);
Culdata_cal_han.city_fa = (data.city_fa);
%%% 经济地位 主观（1-10） 客观（income)
Culdata_cal_han.sostat_subj_ind = (data.subj_sostatus_fam);
Culdata_cal_han.sostat_subj_fam = (data.subj_sostatus_ind);
Culdata_cal_han.income_anu = (data.income_anu);
Culdata_cal_han.income_resp = (data.income_resp);

%%% BIS/BAS
bisbas_reward = (data.bisbas_1) + (data.bisbas_2) + (data.bisbas_3) + (data.bisbas_4);
bisbas_force = (data.bisbas_5) + (data.bisbas_6) + (data.bisbas_7) + (data.bisbas_8);
bisbas_please = (data.bisbas_9) + (data.bisbas_10) + (data.bisbas_11) + (data.bisbas_12) + (data.bisbas_13);
bisbas_bis = (data.bisbas_14) + (data.bisbas_15) + (data.bisbas_16) + (data.bisbas_17) + (data.bisbas_18);
data.bisbas_reward = bisbas_reward;
data.bisbas_force = bisbas_force;
data.bisbas_please = bisbas_please;
data.bisbas_bis = bisbas_bis;
Culdata_cal_han.bisbas_reward = bisbas_reward;
Culdata_cal_han.bisbas_force = bisbas_force;
Culdata_cal_han.bisbas_please = bisbas_please;
Culdata_cal_han.bisbas_bis = bisbas_bis;

%%% hofstede's value
Hofstede_powerdis = 35*(data.Hofstede_7 - data.Hofstede_2) + 25*(data.Hofstede_20 - data.Hofstede_23);
Hofstede_indvcoll = 35*(data.Hofstede_4 - data.Hofstede_1) + 35*(data.Hofstede_9 - data.Hofstede_6);
Hofstede_masvfem = 35*(data.Hofstede_5 - data.Hofstede_3) + 35*(data.Hofstede_8 - data.Hofstede_10);
Hofstede_unceravoi = 40*(data.Hofstede_18 - data.Hofstede_15) + 25*(data.Hofstede_21 - data.Hofstede_24);
Hofstede_lonvsho = 40*(data.Hofstede_13 - data.Hofstede_14) + 25*(data.Hofstede_19 - data.Hofstede_22);
Hofstede_indvres = 35*(data.Hofstede_12 - data.Hofstede_11) + 40*(data.Hofstede_17 - data.Hofstede_16);
data.Hofstede_powerdis = Hofstede_powerdis;
data.Hofstede_indvcoll = Hofstede_indvcoll;
data.Hofstede_masvfem = Hofstede_masvfem;
data.Hofstede_unceravoi = Hofstede_unceravoi;
data.Hofstede_lonvsho = Hofstede_lonvsho;
data.Hofstede_indvres = Hofstede_indvres;

Culdata_cal_han.Hofstede_powerdis = Hofstede_powerdis;
Culdata_cal_han.Hofstede_indvcoll = Hofstede_indvcoll;
Culdata_cal_han.Hofstede_masvfem = Hofstede_masvfem;
Culdata_cal_han.Hofstede_unceravoi = Hofstede_unceravoi;
Culdata_cal_han.Hofstede_lonvsho = Hofstede_lonvsho;
Culdata_cal_han.Hofstede_indvres = Hofstede_indvres;

Culdata_cal_han.tightness_1 = data.tightness_1;
Culdata_cal_han.tightness_2 = data.tightness_2;
Culdata_cal_han.tightness_3 = data.tightness_3;
Culdata_cal_han.tightness_4_r = 5 - data.tightness_4;
Culdata_cal_han.tightness_5 = data.tightness_5;
Culdata_cal_han.tightness_6 = data.tightness_6;

%%
Culdata_cal_han.z_tightness = fzscore(data.tightness_1) + fzscore(data.tightness_2) + fzscore(data.tightness_3) + fzscore(5-data.tightness_4)...
    + fzscore(data.tightness_5) + fzscore(data.tightness_6);
Culdata_cal_han.z_tightness_all = fzscore(data.tightness_1 + data.tightness_2 + data.tightness_3 + data.tightness_4...
    + data.tightness_5 + data.tightness_6);
Culdata_cal_han.z_tightness_all_re = fzscore(data.tightness_1 + data.tightness_2 + data.tightness_3 + 5-data.tightness_4...
    + data.tightness_5 + data.tightness_6);
Culdata_cal_han.z_tightness_all_ex = fzscore(data.tightness_1 + data.tightness_2 + data.tightness_3...
    + data.tightness_5 + data.tightness_6);
%%
Culdata_cal_han.latitude = data.latitude;
Culdata_cal_han.longitude = data.longitude;
province_pool = {'河北省','山西省','辽宁省','吉林省','黑龙江省','江苏省','浙江省','安徽省','福建省','江西省','山东省','河南省','湖北省','湖南省','广东省','海南省','四川省','贵州省','云南省','陕西省','甘肃省','青海省',...
    '内蒙古自治区','广西壮族自治区','西藏自治区','宁夏回族自治区','新疆维吾尔自治区',...
    '北京市','天津市','上海市','重庆市'};
for iii = 1:length(province_pool)
    Culdata_cal_han.provinceID(Culdata_cal_han.province == province_pool(iii))= iii;
end
Culdata_cal_han.provinceID = Culdata_cal_han.provinceID(:);
% cultural control
behavenames = {'Hofstede_indvcoll','Hofstede_indvres','Hofstede_lonvsho','Hofstede_masvfem','Hofstede_powerdis','Hofstede_unceravoi','z_tightness','z_tightness_all_re','z_tightness_all_ex'};
behavenames_tran = {'Hofstede_indvcoll_control','Hofstede_indvres_control','Hofstede_lonvsho_control','Hofstede_masvfem_control','Hofstede_powerdis_control','Hofstede_unceravoi_control','z_tightness_control','z_tightness_all_re_control','z_tightness_all_ex_control'};
Control_X = [ Culdata_cal_han.age(:) Culdata_cal_han.edu_year(:) Culdata_cal_han.edu_fa(:) Culdata_cal_han.edu_ma(:) Culdata_cal_han.sostat_subj_ind(:) Culdata_cal_han.sostat_subj_fam(:) Culdata_cal_han.income_anu(:) Culdata_cal_han.income_resp(:)];
for beid = 1:length(behavenames)
    [~,~,stats] = glmfit(Control_X,Culdata_cal_han.(behavenames{beid}));
    Culdata_cal_han.(behavenames_tran{beid}) = fzscore(stats.resid(:));
end


% 定义北京的经纬度
lat_beijing = 39.9042;
lon_beijing = 116.4074;
china_provinces = {
    'Jilin',            43.8176,  125.3235; % 长春 (Changchun)
    'Sichuan',          30.5728,  104.8801; % 成都 (Chengdu)
    'Ningxia',          38.4664,  106.2782; % 银川 (Yinchuan)
    'Shandong',         36.6512,  117.1201; % 济南 (Jinan)
    'Jiangsu',          32.0603,  118.7965; % 南京 (Nanjing)
    'Hainan',           20.0440,  110.1999; % 海口 (Haikou)
    'Guizhou',          26.6477,  106.6302; % 贵阳 (Guiyang)
    'Qinghai',          36.6170,  101.7781; % 西宁 (Xining)
    'Shanghai',         31.2304,  121.4737; % 上海 (Shanghai)
    'Yunnan',           25.0389,  102.7183; % 昆明 (Kunming)
    'Inner Mongolia',   40.8283,  111.8323; % 呼和浩特 (Hohhot)
    'Beijing',          39.9042,  116.4074; % 北京 (Beijing)
    'Tianjin',          39.3434,  117.3616; % 天津 (Tianjin)
    'Anhui',            31.8206,  117.2272; % 合肥 (Hefei)
    'Shanxi',           37.8706,  112.5489; % 太原 (Taiyuan)
    'Guangdong',        23.1291,  113.2644; % 广州 (Guangzhou)
    'Guangxi',          22.8176,  108.3669; % 南宁 (Nanning)
    'Xinjiang',         43.8256,  87.6168;  % 乌鲁木齐 (Urumqi)
    'Jiangxi',          28.6820,  115.8575; % 南昌 (Nanchang)
    'Hebei',            38.0419,  114.5149; % 石家庄 (Shijiazhuang)
    'Henan',            34.7466,  113.6254; % 郑州 (Zhengzhou)
    'Zhejiang',         30.2741,  120.1551; % 杭州 (Hangzhou)
    'Hubei',            30.5928,  114.3055; % 武汉 (Wuhan)
    'Hunan',            28.2282,  112.9388; % 长沙 (Changsha)
    'Gansu',            36.0611,  103.8343; % 兰州 (Lanzhou)
    'Fujian',           26.0745,  119.2965; % 福州 (Fuzhou)
    'Tibet',            29.6525,  91.1721;  % 拉萨 (Lhasa)
    'Liaoning',         41.8057,  123.4315; % 沈阳 (Shenyang)
    'Chongqing',        29.5630,  106.5516; % 重庆 (Chongqing)
    'Shaanxi',          34.3416,  108.9398; % 西安 (Xi'an)
    'Heilongjiang',     45.7560,  126.6407; % 哈尔滨 (Harbin)
};

% 计算距离函数（哈弗辛公式）
haversine = @(lat1, lon1, lat2, lon2) ...
    6371 * 2 * atan2(sqrt(sin(deg2rad(lat2-lat1)/2).^2 + ...
    cos(deg2rad(lat1)).*cos(deg2rad(lat2)).*sin(deg2rad(lon2-lon1)/2).^2), ...
    sqrt(1-sin(deg2rad(lat2-lat1)/2).^2 - ...
    cos(deg2rad(lat1)).*cos(deg2rad(lat2)).*sin(deg2rad(lon2-lon1)/2).^2));

% 计算每个省会到北京的距离
distances = zeros(size(china_provinces, 1), 1);
for i = 1:size(china_provinces, 1)
    lat = china_provinces{i, 2};
    lon = china_provinces{i, 3};
    distances(i) = haversine(lat_beijing, lon_beijing, lat, lon);
end
distances_fromBJ = distances;

%% Calculate each cultural distance between each city and Beijing
Culnames_used = {'Hofstede_indvcoll_control','Hofstede_indvres_control','Hofstede_lonvsho_control','Hofstede_masvfem_control','Hofstede_powerdis_control','Hofstede_unceravoi_control','z_tightness_all_ex_control'};
p1idx = Culdata_cal_han.province == '北京市';
plist = unique(Culdata_cal_han.province);
variance_reference = Culdata_cal_han.province == plist(end);
Culdis_from_Beijing_all = nan(length(plist),length(Culnames_used));

for ppp = 1:length(plist)
    p2idx = Culdata_cal_han.province == plist(ppp);
    for ccc = 1:length(Culnames_used)
        M_all = mean(Culdata_cal_han.(Culnames_used{ccc})(p1idx|p2idx));
        M_gen1 = mean(Culdata_cal_han.(Culnames_used{ccc})(p1idx));
        M_gen2 = mean(Culdata_cal_han.(Culnames_used{ccc})(p2idx));
        theta_g = (M_gen1 - M_all)^2 + (M_gen2 - M_all)^2;
        theta_t = var(Culdata_cal_han.(Culnames_used{ccc})(variance_reference));
        Culdis_from_Beijing_all(ppp,ccc) = theta_g/theta_t;
    end
end
Culdis_from_Beijing_all_m = mean(Culdis_from_Beijing_all,2);
Culdis_from_Beijing_pc1_m = mean(Culdis_from_Beijing_all(:,[2 6 7]),2);
Culdis_from_Beijing_pc2_m = mean(Culdis_from_Beijing_all(:,[1 3 4 5]),2);

%% bootstrap
Culnames_used = {'Hofstede_indvcoll_control','Hofstede_indvres_control','Hofstede_lonvsho_control','Hofstede_masvfem_control','Hofstede_powerdis_control','Hofstede_unceravoi_control','z_tightness_all_ex_control'};
p1idx = Culdata_cal_han.province == '北京市';
plist = unique(Culdata_cal_han.province);
Culdis_from_Beijing_boot_CI = nan(length(plist),2);
Culdis_from_Beijing_boot_CI_pc1 = nan(length(plist),2);
Culdis_from_Beijing_boot_CI_pc2 = nan(length(plist),2);

Culdis_from_Beijing_boot_sem = nan(length(plist),1);
Culdis_from_Beijing_boot_sem_pc1 = nan(length(plist),1);
Culdis_from_Beijing_boot_sem_pc2 = nan(length(plist),1);
Culdis_from_Beijing_boot_M = nan(length(plist),1);

stimes = draws;
for ppp = 1:length(plist)
    p2idx = Culdata_cal_han.province == plist(ppp);
    group1_pool = find(p1idx);
    group2_pool = find(p2idx);
    culdis_from_Beijing_boot = nan(length(Culnames_used),stimes);
    for ccc=1:length(Culnames_used)
        x=Culdata_cal_han.(Culnames_used{ccc})(group1_pool);
        y=Culdata_cal_han.(Culnames_used{ccc})(group2_pool);
        for first=1:128:stimes
            ii=first:min(first+127,stimes);
            xb=x(randi(numel(x),numel(x),numel(ii)));
            yb=y(randi(numel(y),numel(y),numel(ii)));
            pool=[xb;yb];avg=mean(pool,1);
            culdis_from_Beijing_boot(ccc,ii)=((mean(xb,1)-avg).^2+(mean(yb,1)-avg).^2)./var(pool,0,1);
        end
    end
    Culdis_from_Beijing_boot_M(ppp,:) = mean(mean(culdis_from_Beijing_boot,2));
    Culdis_from_Beijing_boot_CI(ppp,:) = [prctile(mean(culdis_from_Beijing_boot,1),97.5) prctile(mean(culdis_from_Beijing_boot,1),2.5)];
    Culdis_from_Beijing_boot_CI_pc1(ppp,:) = [prctile(mean(culdis_from_Beijing_boot([2 6 7],:),1),97.5) prctile(mean(culdis_from_Beijing_boot([2 6 7],:),1),2.5)];
    Culdis_from_Beijing_boot_CI_pc2(ppp,:) = [prctile(mean(culdis_from_Beijing_boot([1 3 4 5],:),1),97.5) prctile(mean(culdis_from_Beijing_boot([1 3 4 5],:),1),2.5)];
    Culdis_from_Beijing_boot_sem(ppp) = std(mean(culdis_from_Beijing_boot,1));
    Culdis_from_Beijing_boot_sem_pc1(ppp) = std(mean(culdis_from_Beijing_boot([2 6 7],:),1));
    Culdis_from_Beijing_boot_sem_pc2(ppp) = std(mean(culdis_from_Beijing_boot([1 3 4 5],:),1));
end

T=table(string(plist(:)),string(china_provinces(:,1)),distances_fromBJ(:),Culdis_from_Beijing_pc1_m,Culdis_from_Beijing_pc2_m,Culdis_from_Beijing_all_m,Culdis_from_Beijing_boot_sem_pc1,Culdis_from_Beijing_boot_sem_pc2,Culdis_from_Beijing_boot_sem,'VariableNames',{'Province','Plot_label','Distance_from_Beijing_km','CFST_PC1','CFST_PC2','CFST_all','CFST_PC1_bootstrap_SD','CFST_PC2_bootstrap_SD','CFST_all_bootstrap_SD'});
writetable(T,fullfile(out,'EDFig1bc_source_data.csv'),'Encoding','UTF-8');
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','simulation'));
for pc=1:2
    f=figure('Visible','off','Color','w');errorbar(T.Distance_from_Beijing_km,T.(sprintf('CFST_PC%d',pc)),T.(sprintf('CFST_PC%d_bootstrap_SD',pc)),'.','Color',[.88 .56 .51]);xlabel('Distance from Beijing (km)');ylabel('Cultural distance');publication_style(f);
    savefig(f,fullfile(out,sprintf('CFST_PC%d.fig',pc)));close(f);
end
end

