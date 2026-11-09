"""Extract earliest residence and residence-history summaries."""
import os
import pandas as pd
import numpy as np

def extract_primary_residence(input_excel_path):
    print(f"[1/4] 正在读取原始数据文件: {input_excel_path} ...")
    
    if not os.path.exists(input_excel_path):
        raise FileNotFoundError(f"未找到输入文件: {input_excel_path}")
        
    df_geo = pd.read_excel(input_excel_path, sheet_name=0)
    df_mri = pd.read_excel(input_excel_path, sheet_name=1)
    
    print(f"      - 居住史记录行数: {len(df_geo)}")
    print(f"      - 影像被试数量 (ID数): {df_geo['id'].nunique()}")
    print(f"      - MRI实验记录被试数量: {len(df_mri)}")
    
    # 建立 MRI 扫描日期映射
    study_date_map = dict(zip(df_mri['ID'], df_mri['StudyDate']))
    
    print("[2/4] 正在逐个被试提取早期稳定居住地坐标及居住轨迹特征 ...")
    
    results = []
    
    # 按照 ID 分组处理
    for sid, group in df_geo.groupby('id'):
        # 按年份和序号正序排列，保证时间线严谨
        group_sorted = group.sort_values(['year', 'num']).reset_index(drop=True)
        
        years = group_sorted['year'].values
        lons = group_sorted['long_wgs_main address'].values
        lats = group_sorted['lat_wgs_main address'].values
        provs = group_sorted['PROVINCE'].values
        cities = group_sorted['CITY'].values
        
        # 1. 早期稳定居住地（出生/起始年份居住地）
        birth_year = int(years[0])
        early_lon = float(lons[0])
        early_lat = float(lats[0])
        early_prov = str(provs[0])
        early_city = str(cities[0])
        
        # 2. 计算在早期稳定居住地的连续居住年限（从出生到首次迁移前的年数）
        same_loc_mask = (lons == early_lon) & (lats == early_lat)
        diff_indices = np.where(~same_loc_mask)[0]
        
        if len(diff_indices) > 0:
            first_diff_idx = diff_indices[0]
            first_move_year = int(years[first_diff_idx])
            early_stay_years = first_move_year - birth_year
            is_mover = 1
        else:
            first_move_year = np.nan
            early_stay_years = int(years[-1] - birth_year + 1)
            is_mover = 0
            
        # 3. 历年中累计居住时间最长的地点（众数地点）
        loc_tuples = list(zip(lats, lons, provs, cities))
        loc_series = pd.Series(loc_tuples)
        loc_counts = loc_series.value_counts()
        top_loc = loc_counts.index[0]
        top_loc_lat, top_loc_lon, top_loc_prov, top_loc_city = top_loc
        top_loc_years = int(loc_counts.iloc[0])
        
        # 4. 扫描/实验参加年份时的居住地（最后记录年份）
        scan_year = int(years[-1])
        scan_lon = float(lons[-1])
        scan_lat = float(lats[-1])
        scan_prov = str(provs[-1])
        scan_city = str(cities[-1])
        
        results.append({
            'ID': int(sid),
            'StudyDate': int(study_date_map[sid]) if (sid in study_date_map and pd.notna(study_date_map[sid])) else np.nan,
            'BirthYear': birth_year,
            'ScanYear': scan_year,
            'Primary_Longitude_WGS84': early_lon,
            'Primary_Latitude_WGS84': early_lat,
            'Primary_Province': early_prov,
            'Primary_City': early_city,
            'EarlyStayYears': early_stay_years,
            'IsMover': is_mover,
            'FirstMoveYear': first_move_year if pd.notna(first_move_year) else np.nan,
            'TotalRecordedYears': len(group),
            'MostFrequent_Longitude_WGS84': float(top_loc_lon),
            'MostFrequent_Latitude_WGS84': float(top_loc_lat),
            'MostFrequent_Province': str(top_loc_prov),
            'MostFrequent_City': str(top_loc_city),
            'MostFrequent_Years': top_loc_years,
            'Scan_Longitude_WGS84': scan_lon,
            'Scan_Latitude_WGS84': scan_lat,
            'Scan_Province': scan_prov,
            'Scan_City': scan_city
        })
        
    res_df = pd.DataFrame(results)
    
    print(f"[3/4] 提取完成! 共生成 {len(res_df)} 名被试的核心数据。")
    print(f"      - 早期稳定居住平均年限: {res_df['EarlyStayYears'].mean():.2f} 年 (中位数: {res_df['EarlyStayYears'].median():.1f} 年)")
    print(f"      - 早期居住地覆盖省份数: {res_df['Primary_Province'].nunique()} 个")
    print(f"      - 早期居住地覆盖城市数: {res_df['Primary_City'].nunique()} 个")
    
    return res_df
