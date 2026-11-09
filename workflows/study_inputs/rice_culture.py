import pandas as pd
import numpy as np
import scipy.stats as stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Load data
import argparse, json
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('--input', help='Merged province table')
p.add_argument('--province-workbook')
p.add_argument('--indicator-sav')
p.add_argument('--output', required=True, type=Path)
a = p.parse_args()
mapping = {
    '安徽': 'Anhui', '北京': 'Beijing', '福建': 'Fujian', '甘肃': 'Gansu',
    '广东': 'Guangdong', '广西': 'Guangxi', '贵州': 'Guizhou', '海南': 'Hainan',
    '河北': 'Hebei', '河南': 'Henan', '黑龙': 'Heilongjiang', '湖北': 'Hubei',
    '湖南': 'Hunan', '吉林': 'Jilin', '江苏': 'Jiangsu', '江西': 'Jiangxi',
    '辽宁': 'Liaoning', '内蒙': 'Mongolia', '宁夏': 'Ningxia', '青海': 'Qinghai',
    '山东': 'Shandong', '山西': 'Shanxi', '陕西': 'Shaanxi', '上海': 'Shanghai',
    '四川': 'Sichuan', '天津': 'Tianjin', '西藏': 'Tibet', '新疆': 'Xinjiang',
    '云南': 'Yunnan', '浙江': 'Zhejiang', '重庆': 'Chongqing'
}

if a.input and not (a.province_workbook or a.indicator_sav):
    df = pd.read_csv(a.input, encoding='utf-8-sig')
elif a.province_workbook and a.indicator_sav and not a.input:
    provinces = pd.read_excel(a.province_workbook)
    provinces['Location'] = provinces.sname.map(mapping)
    df = provinces.merge(pd.read_spss(a.indicator_sav),on='Location',validate='many_to_one')
else:
    p.error('Supply --input or both --province-workbook and --indicator-sav')

# Keep only province-level variables
cols_to_keep = [
    'growprov', 'PercentPaddy', 'EnvironmentalRiceSuitability', 
    'Hofstede_indvcoll_control', 'Hofstede_indvres_control', 'Hofstede_lonvsho_control', 
    'Hofstede_masvfem_control', 'Hofstede_powerdis_control', 'Hofstede_unceravoi_control', 
    'z_tightness_control'
]
df_prov = df[cols_to_keep].drop_duplicates(subset=['growprov']).dropna(subset=['growprov'])

# Exclude Xinjiang, Tibet, Inner Mongolia
exclude_provs = ['Xinjiang', 'Tibet', 'Inner Mongolia', 'Neimenggu', 'Xizang', '新疆', '西藏', '内蒙古', '内蒙', '内蒙古自治区', '新疆维吾尔自治区', '西藏自治区']
df_prov = df_prov[~df_prov['growprov'].isin(exclude_provs)]
print("Included provinces:", len(df_prov), list(df_prov['growprov']))

# Define cultural variables and rice variables
culture_vars = {
    'Hofstede_indvcoll_control': 'IDV',
    'Hofstede_powerdis_control': 'PDI',
    'Hofstede_masvfem_control': 'MAS',
    'Hofstede_unceravoi_control': 'UAI',
    'Hofstede_lonvsho_control': 'LTO',
    'Hofstede_indvres_control': 'IVR',
    'z_tightness_control': 'TH'
}
rice_vars = {
    'PercentPaddy': 'Percent Paddy',
    'EnvironmentalRiceSuitability': 'Rice Suitability'
}

# 1. Compute correlations for individual cultural dimensions
results = []
for rv, rv_name in rice_vars.items():
    for cv, cv_name in culture_vars.items():
        # drop na
        valid_df = df_prov.dropna(subset=[rv, cv])
        if len(valid_df) < 3: continue
        r, p = stats.pearsonr(valid_df[rv], valid_df[cv])
        results.append({
            'Rice Index': rv_name,
            'Cultural Dimension': cv_name,
            'r': round(r, 3),
            'p-value': round(p, 4),
            'N': len(valid_df)
        })

# 2. PCA on the 7 cultural dimensions
culture_cols = list(culture_vars.keys())
valid_pca_df = df_prov.dropna(subset=culture_cols).copy()

if len(valid_pca_df) >= 3:
    scaler = StandardScaler()
    scaled_culture = scaler.fit_transform(valid_pca_df[culture_cols])
    
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled_culture)
    
    valid_pca_df['Culture_PC1'] = pcs[:, 0]
    valid_pca_df['Culture_PC2'] = pcs[:, 1]
    
    print("\n--- PCA Results (7 Cultural Dimensions) ---")
    print(f"Explained Variance Ratio: PC1 = {pca.explained_variance_ratio_[0]*100:.2f}%, PC2 = {pca.explained_variance_ratio_[1]*100:.2f}% (Total = {np.sum(pca.explained_variance_ratio_)*100:.2f}%)")
    
    loadings_df = pd.DataFrame(
        pca.components_.T, 
        index=[culture_vars[c] for c in culture_cols], 
        columns=['PC1 Loading', 'PC2 Loading']
    )
    print("\nPCA Loadings:")
    print(loadings_df.round(3).to_string())
    
    # Compute correlations for PC1 and PC2
    pca_vars = {
        'Culture_PC1': 'Culture PC1',
        'Culture_PC2': 'Culture PC2'
    }
    
    for rv, rv_name in rice_vars.items():
        for pc_col, pc_name in pca_vars.items():
            valid_df = valid_pca_df.dropna(subset=[rv, pc_col])
            if len(valid_df) < 3: continue
            r, p = stats.pearsonr(valid_df[rv], valid_df[pc_col])
            results.append({
                'Rice Index': rv_name,
                'Cultural Dimension': pc_name,
                'r': round(r, 3),
                'p-value': round(p, 4),
                'N': len(valid_df)
            })

print("\n--- Correlation Results ---")
out=pd.DataFrame(results)
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(dict(columns=list(out.columns),rows=out.values.tolist()),ensure_ascii=False,indent=2),encoding='utf-8')


