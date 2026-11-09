"""Rebuild the provincial psychometric table from individual questionnaire inputs."""
from pathlib import Path
import argparse, json, sys
import numpy as np
import pandas as pd
from scipy import stats
from matlab_source_reader import read_mat

DIMENSIONS={
 'IVR':('Hofstede_indvres',[11,12,16,17]),
 'UAI':('Hofstede_unceravoi',[15,18,21,24]),
 'PDI':('Hofstede_powerdis',[2,7,20,23]),
 'IDV':('Hofstede_indvcoll',[1,4,6,9]),
 'MAS':('Hofstede_masvfem',[3,5,8,10]),
 'LTO':('Hofstede_lonvsho',[13,14,19,22])}
SCORE_TERMS={
 'Hofstede_powerdis':[(35,7,2),(25,20,23)],
 'Hofstede_indvcoll':[(35,4,1),(35,9,6)],
 'Hofstede_masvfem':[(35,5,3),(35,8,10)],
 'Hofstede_unceravoi':[(40,18,15),(25,21,24)],
 'Hofstede_lonvsho':[(40,13,14),(25,19,22)],
 'Hofstede_indvres':[(35,12,11),(40,17,16)]}

def anova_icc(y,groups):
    d=pd.DataFrame({'y':pd.to_numeric(y,errors='coerce'),'group':groups}).dropna()
    if d['group'].dtype==object: d=d[d.group.astype(str).str.len()>0]
    samples=[v.y.to_numpy() for _,v in d.groupby('group')]
    n=len(d); k=len(samples); mu=d.y.mean()
    msb=sum(len(v)*(v.mean()-mu)**2 for v in samples)/(k-1)
    msw=sum(((v-v.mean())**2).sum() for v in samples)/(n-k)
    F,p=stats.f_oneway(*samples)
    icc1=(msb-msw)/(msb+(n/k-1)*msw)
    return {'n':n,'groups':k,'F':float(F),'p':float(p),'ICC1':float(icc1),'ICC2':float((msb-msw)/msb)}

def rwg(data,items,maximum):
    null=(maximum**2-1)/12
    rows=[]
    for name,g in data.dropna(subset=['from_province']).groupby('from_province'):
        if len(g)<2: continue
        mv=g[items].apply(pd.to_numeric,errors='coerce').var(ddof=1).mean()
        j=len(items); v=0. if mv>null else j*(1-mv/null)/(j*(1-mv/null)+mv/null)
        rows.append({'province':str(name),'n':len(g),'item_variance':float(mv),'rwg':float(v)})
    return float(np.mean([v['rwg'] for v in rows])),rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--questionnaire',required=True);p.add_argument('--culture-mat',required=True);p.add_argument('--tl-null-options',type=int,choices=[5,6],required=True);p.add_argument('--output-dir',required=True);a=p.parse_args()
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    df=pd.read_excel(a.questionnaire)
    score_checks=[]
    for name,terms in SCORE_TERMS.items():
        score=sum(weight*(pd.to_numeric(df[f'Hofstede_{i}'])-pd.to_numeric(df[f'Hofstede_{j}'])) for weight,i,j in terms)
        saved=pd.to_numeric(df[name]);err=float(np.nanmax(np.abs(score-saved)))
        np.testing.assert_allclose(score,saved,rtol=1e-12,atol=1e-12,equal_nan=True)
        score_checks.append({'dimension_field':name,'n':len(score),'max_error':err})
        df[name]=score
    for i in range(1,7): df[f'tightness_{i}']=pd.to_numeric(df[f'tightness_{i}'],errors='coerce')
    rows=[];group_rows=[];variants=[]
    for name,(var,items) in DIMENSIONS.items():
        f=anova_icc(df[var],df.from_province)
        ic=anova_icc(df[var],df.province)
        mean,parts=rwg(df,[f'Hofstede_{i}' for i in items],5)
        rows.append({'dimension':name,'alpha':np.nan,'rwg_mean':mean,'ICC2':ic['ICC2'],'F':f['F'],'p':f['p'],'anova_n':f['n'],'icc_n':ic['n'],'anova_groups':f['groups'],'icc_groups':ic['groups']})
        group_rows.extend([{'dimension':name,**v} for v in parts])
    for items in [[1,2,3,5,6],[1,2,3,4,5,6]]:
        y=df[[f'tightness_{i}' for i in items]].mean(axis=1)
        f=anova_icc(y,df.from_province);ic=anova_icc(y,df.province)
        mean,parts=rwg(df,[f'tightness_{i}' for i in items],6)
        variants.append({'items':items,'anova':f,'province_icc':ic,'rwg_mean':mean})
        if len(items)==6: group_rows.extend([{'dimension':'TL',**v} for v in parts])
    c=read_mat(a.culture_mat,'Culdata')
    x=np.column_stack([np.asarray(c[f'tightness_{i}']).ravel() for i in [1,2,3,5,6]])
    total=np.var(x.sum(axis=1),ddof=1); alpha=5/4*(total-np.var(x,axis=0,ddof=1).sum())/total
    alpha_mat=float(alpha)
    x=df[[f'tightness_{i}' for i in [1,2,3,5,6]]].to_numpy()
    total=np.var(x.sum(axis=1),ddof=1);alpha=5/4*(total-np.var(x,axis=0,ddof=1).sum())/total
    f=variants[0]['anova'];ic=variants[0]['province_icc']
    tl_rwg,parts=rwg(df,[f'tightness_{i}' for i in [1,2,3,5,6]],a.tl_null_options)
    group_rows=[r for r in group_rows if r['dimension']!='TL']+[{'dimension':'TL',**v} for v in parts]
    rows.append({'dimension':'TL','alpha':float(alpha),'rwg_mean':tl_rwg,'ICC2':ic['ICC2'],'F':f['F'],'p':f['p'],'anova_n':f['n'],'icc_n':ic['n'],'anova_groups':f['groups'],'icc_groups':ic['groups']})
    pd.DataFrame(rows).to_csv(out/'EDFig1d_statistics.csv',index=False)
    pd.DataFrame(group_rows).to_csv(out/'EDFig1d_province_agreement.csv',index=False)
    settings={'questionnaire':str(Path(a.questionnaire).resolve()),'culture_mat':str(Path(a.culture_mat).resolve()),'icc_group_column':'province','anova_rwg_group_column':'from_province','TL_icc_items':[1,2,3,5,6],'TL_anova_rwg_items':[1,2,3,5,6],'TL_null_options':a.tl_null_options,'TL_observed_response_min':float(np.min(x)),'TL_observed_response_max':float(np.max(x)),'alpha_items':[1,2,3,5,6],'alpha_input':'questionnaire','alpha_n':len(x),'alpha_missing_cells':int(np.isnan(x).sum())}
    (out/'psychometric_results.json').write_text(json.dumps({'statistics':rows,'tightness_variants':variants,'culture_mat_alpha':alpha_mat,'item_score_verification':score_checks,'configuration':settings},indent=2,allow_nan=True),encoding='utf8')
    print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__': main()
