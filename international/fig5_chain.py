"""Raw cortex, within-unit ComBat, endpoint models and Fig. 5 comparisons."""
from pathlib import Path
import argparse,json,os,sys
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
import numpy as np
import pandas as pd
from scipy.io import savemat
from run_example import z,unit,make,artifact_endpoint,fit,design
from endpoint_permutations import permutation_test
from combat import combat,combat_plan
from model_spec import NATIONS_BLOCKS,USA_MACROS
from wald import wald_statistics
from webb_intervals import webb_intervals
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'figures'))
from render_maps import render

MODELS={'PC1':('Nations','eigenmode2','China_PC1',-1),'PC2':('Nations','SA','China_PC2',1),'USA_TL':('USA','eigenmode2','China_TL',1)}

def align(frame,key,names):
    if frame[key].isna().any() or frame[key].duplicated().any():raise ValueError('Missing or duplicated '+key)
    missing=set(names)-set(frame[key]);
    if missing:raise ValueError(f'Missing {key}: {sorted(missing)}')
    return frame.set_index(key).loc[list(names)]

def read_cortex(folder,cohort):
    folder=folder/cohort
    meta=pd.read_csv(folder/'participants.csv',dtype={'row_key':str,'unit':str,'batch':str})
    if meta[['row_key','unit','batch']].isna().any().any() or meta.row_key.duplicated().any():raise ValueError('Participant keys must be complete and unique')
    table=align(pd.read_csv(folder/'cortex.csv',dtype={'row_key':str}),'row_key',meta.row_key)
    raw=table.to_numpy(float);g,names=pd.factorize(meta.unit,sort=False);batch,_=pd.factorize(meta.batch,sort=False)
    d=dict(g=g,group_names=np.array(names,str),roi_names=np.array(table.columns,str),age=meta.age.to_numpy(float),female=meta.female.to_numpy(float),tiv=meta.tiv.to_numpy(float),raw=raw)
    for key in ['age_observed','sex_observed']:
        flags=meta[key].to_numpy(float)
        if not np.isin(flags,[0,1]).all():raise ValueError('Expected 0/1 '+key)
        d[key]=flags.astype(bool)
    if not all(np.isfinite(d[k]).all() for k in ['age','female','tiv','raw']):raise ValueError('Nonfinite individual inputs')
    d.update(sex_female_design=d['female'],sex_other_design=np.zeros(len(g)))
    d['harm']=combat(raw,combat_plan(None,d,g,batch));return d

def model_inputs(raw,folder,model):
    cohort,axis,china,_=MODELS[model];names=raw['group_names'];roi=raw['roi_names'];d=dict(raw)
    units=align(pd.read_csv(folder/cohort/'units.csv',dtype={'unit':str}),'unit',names)
    macro=sum(NATIONS_BLOCKS.values(),[]) if cohort=='Nations' else USA_MACROS
    d['macro']=units[macro].to_numpy(float);c=units[model].to_numpy(float);d['C']=z(c)
    targets=align(pd.read_csv(folder/'targets.csv'),'roi',roi);d['targets']=targets[[axis,china]].to_numpy(float)
    with np.load(folder/cohort/'pool.npz') as f:
        if not np.array_equal(f['unit_names'],names):raise ValueError('Artifact pool rows use a different unit order')
        d['pool']=f['pool']
    if d['pool'].dtype.kind not in 'iu' or d['pool'].ndim!=2 or not np.all(np.sort(d['pool'],axis=1)==np.arange(len(names))):raise ValueError('Invalid cultural-unit permutations')
    d['X'],labels=design(d,cohort=='Nations');d['Y']=z(d['harm']);return d,c,units,labels

def save_mat_npz(folder,stem,data):
    folder.mkdir(parents=True,exist_ok=True);np.savez_compressed(folder/(stem+'.npz'),**data)
    savemat(folder/(stem+'.mat'),{k:(v.astype(object) if isinstance(v,np.ndarray) and v.dtype.kind=='U' else v) for k,v in data.items()},do_compression=True)

def comparisons(folder,model,roi,beta,targets,spin,null_r,webb_r,null_method):
    if spin.ndim!=2 or spin.shape[0]!=len(roi) or spin.dtype.kind not in 'iu' or np.any((spin<0)|(spin>=len(roi))):raise ValueError('Invalid zero-based ROI x rotation indices')
    if targets.shape!=(len(roi),2) or null_r.ndim!=2 or webb_r.ndim!=2 or null_r.shape[1]!=2 or webb_r.shape[1]!=2 or not len(null_r) or len(webb_r)<2:raise ValueError('Two aligned comparisons required')
    for j,name in enumerate(MODELS[model][1:3]):
        x=targets[:,j];y=beta;r=float(np.corrcoef(x,y)[0,1]);sn=unit(x[spin]).T@unit(y[:,None]);sn=sn.ravel();pn=null_r[:,j];wb=webb_r[:,j]
        d=dict(roi=roi,x=x,y=y,r=r,spin_null=sn,permutation_null=pn,webb_r=wb,CI95=np.percentile(wb,[2.5,97.5]),p_spin=(abs(sn)>=abs(r)).mean(),p_permutation=(abs(pn)>=abs(r)).mean(),target_name=name,null_method=null_method)
        save_mat_npz(folder,f'Fig5_{model}_comparison{j+1}',d)

def geography(folder,model,units,c,score,abbreviations):
    cohort,_,_,sign=MODELS[model];display=c*sign
    if np.corrcoef(score,c)[0,1]<1-1e-10:raise ValueError('Culture-map values and fitted score differ in orientation or unit ordering')
    save_mat_npz(folder,f'Fig5_{model}_geography',dict(names=units,abbreviations=abbreviations,culture_score=score,raw_score=c,display_sign=sign,display_score=display))
    regions='world' if cohort=='Nations' else 'usa';assets=json.loads((ROOT/f'assets/geography/{regions}.json').read_text(encoding='utf-8'));known={f['properties']['name'] for f in assets['features']}
    aliases={'SouthKorea':['South Korea','Republic of Korea'],'Czechia':['Czech Republic'],'Russia':['Russian Federation'],'United States':['United States of America']}
    keys=[]
    for u in units:
        key=next((v for v in [u.replace('_',' '),*aliases.get(u,[])] if v in known),None)
        if key is None:raise ValueError('Missing geographic boundary: '+u)
        keys.append(key)
    panel=dict(region=regions,title={'PC1':'Culture PC1','PC2':'Culture PC2','USA_TL':'Tightness–looseness'}[model],region_values=dict(zip(keys,np.round(display,2))),color_limits=[-4,4] if cohort=='Nations' else [20,70],colors=['#ffffff','#9bbb59','#52642e'] if cohort=='Nations' else ['#edf1f9','#8da4c9','#4163a0'])
    if cohort=='USA':
        panel['abbreviations']=dict(zip(keys,abbreviations));panel['exclude_regions']=['Alaska','Hawaii','Puerto Rico']
        panel['insets']=[dict(position=[x,0,.15,.15],panel=dict(region='usa',only_regions=[state],title=label,
            region_values=panel['region_values'],colors=panel['colors'],color_limits=panel['color_limits'],colorbar=False)) for x,state,label in [(0,'Alaska','AK'),(.15,'Hawaii','HI')]]
    render(dict(panels=[panel]),folder/f'Fig5_{model}_map')

def run(folder,output,permutations,webb_draws,seed=42,harmonized=None,models=None,resume=False,permutation_indices=None):
    if permutations<2 or webb_draws<2:raise ValueError('At least two permutations and Webb draws required')
    if output.exists() and not resume:raise FileExistsError('Choose a new run directory or use --resume with unchanged inputs')
    models=list(MODELS) if models is None else models
    supplied=None
    if permutation_indices is not None:
        with np.load(permutation_indices) as f:supplied=dict(f)
    for cohort in ['Nations','USA']:
        if not any(MODELS[m][0]==cohort for m in models):continue
        if harmonized:
            with np.load(harmonized/(cohort+'_harmonized.npz')) as f:raw=dict(f)
        else:raw=read_cortex(folder,cohort)
        local=output/'local';local.mkdir(parents=True,exist_ok=True);np.savez_compressed(local/(cohort+'_harmonized.npz'),**raw)
        for model,(group,_,_,_) in MODELS.items():
            if group!=cohort or model not in models:continue
            d,c,units,labels=model_inputs(raw,folder,model);roi=d['roi_names'];names=d['group_names'];src=output/'source'
            with np.load(folder/'spins.npz') as f:
                if not np.array_equal(roi,f['roi_names']):raise ValueError('Spatial rotation ROI order differs')
                spins=f['indices']
            state=make(d['X'],d['Y'],d);a,_=artifact_endpoint(state,d['C'],d['pool'],'screen_r001');beta,lam,clean=fit(a,d['C'],a['k'],True);assert a['k']==len(names)-1
            se,t=wald_statistics(d,beta,lam,clean);obs=(unit(beta[:,None]).T@unit(d['targets']))[0];j=list(MODELS).index(model)
            orders=None
            if supplied is not None:
                if model not in supplied or model+'_unit_names' not in supplied:raise ValueError('Missing permutation indices or unit names for '+model)
                if not np.array_equal(supplied[model+'_unit_names'],names):raise ValueError('Permutation unit order differs for '+model)
                orders=supplied[model]
            null=permutation_test(state,d,obs,permutations,seed+j,local/(model+'_permutations'),orders=orders);webb=webb_intervals(d,webb_draws,2026102201+j);wr=unit(webb['map_draws'].T).T@unit(d['targets'])
            save_mat_npz(src,f'Fig5_{model}_brain',dict(roi=roi,beta=beta,se=se,t=t,endpoint_k=a['k']))
            comparisons(src,model,roi,beta,d['targets'],spins,null['null_r'],wr,null['null_method'])
            abbr=units['abbreviation'].to_numpy(str) if cohort=='USA' else names;geography(src,model,names,c,d['C'],abbr)
            np.savez_compressed(local/(model+'_fit.npz'),X=d['X'],labels=np.array(labels),clean=clean,lambda_ratio=lam,**null,webb_weights=webb['weights'],webb_maps=webb['map_draws'])
            print(model,'endpoint',a['k'],'r',obs.tolist(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--harmonized',type=Path);p.add_argument('--permutations',type=int,default=10000);p.add_argument('--webb-draws',type=int,default=1000);p.add_argument('--seed',type=int,default=42)
    p.add_argument('--models',choices=MODELS,nargs='+',default=list(MODELS));p.add_argument('--resume',action='store_true')
    p.add_argument('--permutation-indices',type=Path,help='NPZ with model-specific draw-by-unit indices and unit names')
    args=p.parse_args();run(args.inputs,args.output,args.permutations,args.webb_draws,args.seed,args.harmonized,args.models,args.resume,args.permutation_indices)
