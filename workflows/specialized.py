"""Figure-specific calculation and source-data adapters."""
import itertools
import numpy as np
from scipy.stats import rankdata,norm,pearsonr,t
from panel_data import panel,labels
from map_statistics import standardize,commonality,spin_test,bh_fdr
from cultural_pca import cultural_pca
from group_summaries import summaries
from special_statistics import traveling
from psychometrics import anova_icc,rwg
from cultural_pca import cronbach_alpha

def heat(figure,name,matrix,row,col=None,p=None):
    a=np.asarray(matrix);low=float(np.nanmin(a));high=float(np.nanmax(a))
    if low==high:low-=.5;high+=.5
    d=panel(figure,name,'heatmap',matrix=a,labels=np.asarray(row,object),xlabels=np.asarray(row if col is None else col,object),clim=[low,high])
    if p is not None:d['p']=p
    return d

def province_summary(d):
    x=summaries(d['x'],d['province']);y=summaries(d['y'],d['province']);P=[]
    if np.any(x['n']<2):raise ValueError('SEM requires at least two rows per province')
    for j in range(x['mean'].shape[1]):
        r,p=pearsonr(x['mean'][:,j],y['mean'][:,j])
        P.append(panel('Figure2',f'Province_PC{j+1}','error_scatter',x=x['mean'][:,j],y=y['mean'][:,j],xse=x['sem'][:,j],se=y['sem'][:,j],labels=x['labels'],color=[.72,.5,.46],stats=f'r = {r:.3f}, p = {p:.3f}',xlabel='Culture',ylabel='Adjusted cortical measure'))
    return dict(province=x['labels'],n=x['n'],culture_mean=x['mean'],culture_sem=x['sem'],brain_mean=y['mean'],brain_sem=y['sem']),P,None

def demographics(d):
    values=np.asarray(d['continuous']);groups=np.asarray(d['group']).astype(str);names=np.unique(groups);rows=[];labs=[]
    for group in names:
        for j,label in enumerate(d['continuous_labels']):
            x=values[groups==group,j];q=np.quantile(x,[.25,.5,.75]);rows.append([len(x),x.mean(),x.std(ddof=1),*q]);labs.append(group+': '+str(label))
    a=np.asarray(rows);P=[panel('ExtendedDataFigure1','Demographics','table',matrix=a,labels=np.array(labs,object),columns=np.array(['N','Mean','SD','Q25','Median','Q75'],object))]
    out=dict(continuous_summary=a,row_labels=np.array(labs))
    if 'categorical' in d:
        rows=[];labs=[]
        for group in names:
            vals,counts=np.unique(np.asarray(d['categorical'])[groups==group].astype(str),return_counts=True)
            for v,n in zip(vals,counts):rows.append([n,100*n/(groups==group).sum()]);labs.append(group+': '+v)
        out.update(categorical_summary=np.array(rows),category_labels=np.array(labs));P.append(panel('ExtendedDataFigure1','Demographic_categories','table',matrix=np.array(rows),labels=np.array(labs,object),columns=np.array(['N','Percent'],object)))
    return out,P,None

def ecological_comparison(d):
    source,P,_=province_summary(d)
    for p in P:p.update(folder='ExtendedDataFigure1',xlabel=str(d.get('x_label','Culture measure')),ylabel=str(d.get('y_label','Reference measure')))
    return source,P,None

def psychometric_table(d):
    import pandas as pd
    scores=np.asarray(d['dimension_scores']);items=np.asarray(d['items']);mask=np.asarray(d['item_mask'],bool);rows=[]
    for j in range(scores.shape[1]):
        f=anova_icc(scores[:,j],d['anova_group']);ic=anova_icc(scores[:,j],d['icc_group']);X=items[:,mask[j]]
        frame=pd.DataFrame(X);frame['from_province']=d['anova_group'];agreement,_=rwg(frame,list(range(X.shape[1])),float(d['response_max'][j]))
        alpha=cronbach_alpha(X) if d['compute_alpha'][j] else np.nan
        rows.append([alpha,agreement,ic['ICC2'],f['F'],f['p'],f['n']])
    matrix=np.array(rows);columns=np.array(['Alpha','rwg','ICC2','F','p','N'],object)
    return dict(statistics=matrix,labels=d['dimension_labels'],columns=columns),[panel('ExtendedDataFigure1','Psychometrics','table',matrix=matrix,labels=d['dimension_labels'],columns=columns)],None

def brain_icc(d):
    Y=np.asarray(d['adjusted_brain']);groups=d['province'];one=[];three=[]
    for a in Y:
        s=[anova_icc(v,groups) for v in a.T];one.append([v['ICC1'] for v in s]);three.append([v['ICC2'] for v in s])
    out=dict(icc1=np.asarray(one),icc3=np.asarray(three),measure_labels=d['measure_labels']);P=[]
    for key in ['icc1','icc3']:
        values=out[key].mean(1);P.append(panel('ExtendedDataFigure1','Brain_'+key,'bars',y=values,labels=d['measure_labels'],color=np.tile([.5,.65,.75],(len(Y),1)),ylabel='Mean '+key.upper()))
    return out,P,None

def gradients(d):
    d=dict(d)
    if 'mode2_display' not in d:d['mode2_display']=-np.asarray(d['DK308_modes_asym'])[:,1]
    if 'fcgradient1' not in d:d['fcgradient1']=np.asarray(d['fcgradientsM'])[:,0]
    fields=['SAaxis','mode2_display','xu_evo','reardon_dev','intervar','cogPC1','myelin','gene_PC1','fcgradient1']
    a=np.column_stack([d[k] for k in fields]);names=np.asarray(d.get('gradient_labels',fields),object);c=np.corrcoef(a.T)
    return dict(maps=a,names=names,correlations=c),[heat('ExtendedDataFigure3','Gradient_correlations',c,names)],None

def commonality_panels(d):
    target=np.atleast_2d(d['mover_maps'])-np.atleast_2d(d['stayer_maps']);P=[];out={}
    for j,y in enumerate(target):
        parts,b,pred=commonality(y,d['deculturation_maps'][j],d['acculturation_maps'][j],d['relocation_maps'][j]);prefix=f'pc{j+1}_'
        out.update({prefix+'target':y,prefix+'fitted':pred,prefix+'bar_values':parts,prefix+'coefficients':b})
        P.append(panel('Figure3',f'Commonality_PC{j+1}','bars',y=parts[[0,2,1,3,4]],labels=np.array(['C','Shared','A','R','Total'],object),color=np.tile([.55,.65,.8],(5,1)),ylabel='Explained variance'))
        P.append(panel('Figure3',f'Commonality_fit_PC{j+1}','identity_scatter',x=y,y=pred,xlabel='Relocator minus stayer',ylabel='Fitted difference'))
    return out,P,None

def best_order(effects):
    # Rows are ROIs, columns are tasks; signed lower-triangle RDMs.
    E=np.asarray(effects,float);R,T=E.shape
    if T>8:raise ValueError('Exhaustive ordering is limited to eight tasks')
    i,j=np.tril_indices(T,-1);v=E[:,i]-E[:,j]
    orders=np.array(list(itertools.permutations(range(T))))
    ranks=np.argsort(orders,axis=1)+1;model=ranks[:,i]-ranks[:,j]
    rv=standardize(rankdata(v,axis=1).T).T;rm=standardize(rankdata(model,axis=1).T).T
    rho=rv@rm.T/(len(i)-1);scores=abs(rho).mean(0);k=np.flatnonzero(scores>=scores.max()-1e-12)[0]
    x=ranks[k]-ranks[k].mean();slopes=E@x/(x@x)
    return dict(best_order=orders[k],rsa_map=rho[:,k],slope=slopes,all_order_rho=rho,all_fit_scores=scores,orders=orders,
                ordered_task_effect=E[:,orders[k]],model_rdm=ranks[k,:,None]-ranks[k,None,:])

def hierarchy(d):
    branches=np.asarray(d['task_effects']);names=np.asarray(d['branch_labels']);tasks=np.asarray(d['task_labels']);P=[];out={}
    if branches.ndim!=3 or len(names)!=len(branches):raise ValueError('task_effects must be branch x ROI x task')
    for j,effects in enumerate(branches):
        fit=best_order(effects);prefix=f'branch{j+1}_';out.update({prefix+k:v for k,v in fit.items()});order=fit['best_order']
        for r in np.asarray(d['example_roi'],int):
            y=fit['ordered_task_effect'][r];rdm=y[:,None]-y[None,:]
            P.append(heat('Figure4',f'Branch{j+1}_ROI{r}_RDM',rdm,tasks[order]))
            P.append(panel('ExtendedDataFigure12',f'Branch{j+1}_ROI{r}_slope','task_scatter',x=np.arange(1,len(tasks)+1),y=y,labels=tasks[order],color=[.6,.45,.65]))
        r,null,p=spin_test(fit['rsa_map'],fit['slope'],d['spin_indices']);out.update({prefix+'rsa_slope_r':r,prefix+'rsa_slope_null':null,prefix+'rsa_slope_p':p})
        P.append(panel('ExtendedDataFigure12',f'Branch{j+1}_RSA_slope','scatter',x=standardize(fit['rsa_map']),y=standardize(fit['slope']),lh=np.asarray(d['left_hemisphere'],bool),xlabel='RSA',ylabel='Task slope',stat_null=null,stat_observed=r,stat_label='r',stats=f'r = {r:.3f}, p = {p:.3f}'))
    out.update(branch_labels=names,task_labels=tasks)
    return out,P,None

def correlation_table(d,figure='ExtendedDataFigure8'):
    x=np.asarray(d['group_values']);n=x.shape[1];c=np.eye(n);p=np.zeros((n,n))
    for i in range(n):
        for j in range(i):c[i,j],p[i,j]=pearsonr(x[:,i],x[:,j]);c[j,i]=c[i,j];p[j,i]=p[i,j]
    return dict(correlations=c,p=p,labels=d['labels']),[heat(figure,figure+'_correlations',c,d['labels'],p=p)],None

def effect_pca(d):
    fit=cultural_pca(d['dimension_effects']);scores=fit['scores'];reference=d['culture_effects'];P=[];out=dict(fit)
    matches=np.asarray(d['reference_component'],int)
    for j,k in enumerate(matches):
        sign=1 if np.corrcoef(scores[:,k],reference[:,j])[0,1]>=0 else -1
        # Orient each component once using its first matched culture map.
        if k not in matches[:j]:scores[:,k]*=sign;fit['loadings'][:,k]*=sign
        r,null,p=spin_test(scores[:,k],reference[:,j],d['spin_indices']);out[f'model{j+1}_r']=r;out[f'model{j+1}_null']=null;out[f'model{j+1}_p']=p
        P.append(panel('ExtendedDataFigure15',f'Effect_PCA_PC{j+1}','scatter',x=standardize(reference[:,j]),y=standardize(scores[:,k]),lh=np.asarray(d['left_hemisphere'],bool),xlabel='Cultural effect',ylabel=f'Effect PC{k+1}',stat_null=null,stat_observed=r,stat_label='r',stats=f'r = {r:.3f}, p = {p:.3f}'))
    P.append(heat('ExtendedDataFigure15','Effect_PCA_loadings',fit['loadings'],d['dimension_labels'],labels('PC',fit['loadings'].shape[1])))
    out['scores']=scores;return out,P,None

def network(d,repetitions,seed):
    groups=['HH','LL','HL','LH','MHH','MLL'];a=[np.asarray(d[g])[:,0]-np.asarray(d[g])[:,1] for g in groups];n=min(len(a[2]),len(a[3]));rng=np.random.default_rng(seed)
    H=np.array([rng.choice(a[0],n,replace=True).mean() for _ in range(repetitions)]);L=np.array([rng.choice(a[1],n,replace=True).mean() for _ in range(repetitions)]);null=H-L
    cross=a[2].mean()-a[3].mean();same=a[4].mean()-a[5].mean()
    # Original network panel uses the smaller strict tail with a 1/B floor.
    tail=lambda v:max(1/repetitions,min(np.mean(null>v),np.mean(null<v)))
    out=dict(null_H=H,null_L=L,null_difference=null,cross=cross,same=same,p_cross=tail(cross),p_same=tail(same),means=np.array([v.mean() for v in a]),group_n=np.array([len(v) for v in a]),groups=np.array(groups,object),sample_size=n)
    return out,[panel('ExtendedDataFigure11','Network_bootstrap','paired_histogram',**out)],None

def resampling(d):
    full=np.asarray(d['full_maps']);draws=np.asarray(d['refit_maps']);P=[];out={}
    # series x draws x proportions x ROI
    if draws.ndim!=4 or draws.shape[0]!=len(full):raise ValueError('refit_maps must be series x draws x proportions x ROI')
    r=np.empty(draws.shape[:3])
    for j in range(len(full)):
        ref=standardize(full[j]);r[j]=standardize(draws[j].reshape(-1,len(ref)).T).T.reshape(draws[j].shape)@ref/(len(ref)-1)
    series=np.stack([np.median(r,axis=1),np.quantile(r,.025,axis=1),np.quantile(r,.975,axis=1)],axis=-1)
    out=dict(correlation_draws=r,series=series,proportions=d['proportions'],series_labels=d['series_labels'])
    P.append(panel('ExtendedDataFigure9','Refit_recovery','ribbons',x=d['proportions'],series=series,labels=d['series_labels']))
    return out,P,None

def word_layout(d):
    counts=np.asarray(d['word_counts'],float);raw=np.asarray(d['scores'],float);y=np.asarray(d['validity_score']);N=len(y)
    if counts.shape!=(N,20) or raw.shape!=(N,3):raise ValueError('Expected 20 word counts and difference/tightness/looseness scores')
    np.testing.assert_allclose(raw[:,0],raw[:,1]-raw[:,2],atol=1e-12)
    if np.any(counts<0) or np.any(counts.sum(1)<=0):raise ValueError('Nonnegative word counts with positive totals required')
    variance=[];loadings=[];pc=[]
    for j in range(2):
        x=counts[:,j*10:(j+1)*10]/counts.sum(1)[:,None];x-=x.mean(0);u,s,v=np.linalg.svd(x,full_matrices=False);score=u[:,0]*s[0];sign=1 if score@raw[:,j+1]>0 else -1
        v[0]*=sign;pc.append(score*sign);variance.append(100*s*s/(s*s).sum());loadings.append(v.T)
    x=np.column_stack([norm.ppf((rankdata(v,method='dense')-.375)/(N+.25)) for v in raw.T]);stats=np.array([pearsonr(v,y) for v in x.T]);fits=[]
    for v in x.T:
        slope,intercept=np.polyfit(v,y,1);xx=np.linspace(v.min()-.05*np.ptp(v),v.max()+.05*np.ptp(v),100);yy=slope*xx+intercept
        se=np.sqrt(np.sum((y-(slope*v+intercept))**2)/(N-2));delta=t.ppf(.975,N-2)*se*np.sqrt(1/N+(xx-v.mean())**2/np.sum((v-v.mean())**2));fits.append(np.c_[xx,yy,yy-delta,yy+delta])
    out=dict(word_counts=counts,word_labels=np.asarray(d['word_labels'],object),dictionary=np.asarray(d['dictionary'],object),province_names=np.asarray(d['province_names'],object),scores=raw,x=x,y=y,r=stats[:,0],p=stats[:,1],explained=np.array(variance),loadings=np.array(loadings),word_frequency_pca_scores=np.array(pc).T,fit=np.array(fits))
    return out,[],'reproduce_suppfig1'

def traveling_unshared(d):
    import re
    a=np.asarray(d['participant1']);b=np.asarray(d['participant2']);na=np.asarray(d['centers1']).astype(str);nb=np.asarray(d['centers2']).astype(str)
    if len(set(na))!=len(na) or len(set(nb))!=len(nb):raise ValueError('Center identifiers must be unique per participant')
    key=lambda s:[int(t) if t.isdigit() else t for t in re.split(r'(\d+)',s)]
    common=np.array(sorted(set(na)&set(nb),key=key));all_names=np.r_[common,sorted((set(na)|set(nb))-set(common),key=key)]
    aa=a[[list(na).index(n) for n in common]];bb=b[[list(nb).index(n) for n in common]];out=traveling(aa,bb)
    dev=[100*(v/np.median(v,axis=0)-1) for v in [a,b]];maps=np.full((2,len(all_names),a.shape[1]),np.nan)
    for j,(v,names) in enumerate(zip(dev,[na,nb])):
        for k,n in enumerate(names):maps[j,list(all_names).index(n)]=v[k]
    spatial=np.array([np.corrcoef(maps[0,k],maps[1,k])[0,1] for k in range(len(common))]);out.update(heatmaps=maps,spatial_r=spatial,median_r=np.median(spatial),centers=all_names.astype(object),common_centers=common.astype(object),scan_counts=[len(a),len(b)],span_medians=[np.median(np.ptp(v,axis=1)) for v in dev])
    return out,[],'reproduce_suppfig2'
