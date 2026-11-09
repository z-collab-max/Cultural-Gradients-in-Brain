"""Assemble numeric figure inputs into source data and MATLAB panel recipes."""
from pathlib import Path
import argparse,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analyses'))
from panel_data import panel,save,labels
from map_panels import assemble_maps
from special_statistics import cultural_rsa,relocation_bootstrap,traveling
from cultural_pca import leave_one_dimension_out,joint_pca_axis_angles
from psychometrics import anova_icc

MAPS=['Figure2','Figure5','ExtendedDataFigure4','ExtendedDataFigure5','ExtendedDataFigure7','ExtendedDataFigure8','ExtendedDataFigure14']
FIGURES=['Figure'+str(i) for i in range(1,6)]+['ExtendedDataFigure'+str(i) for i in range(1,16)]+['SupplementaryFigure'+str(i) for i in range(1,5)]

def assemble(figure,d,repetitions=10000,seed=42):
    for key,value in d.items():
        a=np.asarray(value)
        if a.dtype.kind in 'fc' and not np.isfinite(a).all():raise ValueError('Nonfinite input: '+key)
    panels=[];source=dict(d)
    def add(name,kind,**kw):panels.append(panel(figure,figure+'_'+name,kind,**kw))
    def heat(name,a,row=None,col=None):
        a=np.asarray(a);lo=float(np.min(a));hi=float(np.max(a))
        if lo==hi:lo-=.5;hi+=.5
        add(name,'heatmap',matrix=a,labels=labels('Row',len(a)) if row is None else row,
            xlabels=labels('Column',a.shape[1]) if col is None else col,clim=np.array([lo,hi]))
    if figure in MAPS:
        return assemble_maps(figure,d,repetitions,seed,with_ci=figure in ['Figure2','Figure5','ExtendedDataFigure14'])
    if figure=='ExtendedDataFigure6':
        computed=cultural_rsa(d['cortex'],d['culture'],d['province'],repetitions,seed,str(d['correlation_method']))
        mapped,panels=assemble_maps(figure,{**d,**computed},repetitions,seed)
        return {**computed,**mapped},panels
    if figure=='ExtendedDataFigure3':
        from specialized import gradients
        source,panels,_=gradients(d);return source,panels
    if figure=='Figure1':
        ids=d['province'] if figure=='Figure1' else d['center']
        names,counts=np.unique(ids,return_counts=True);source=dict(group=names,count=counts)
        add('counts','bars',y=counts,labels=names,color=np.tile([.45,.58,.72],(len(counts),1)),ylabel='Participants')
        if figure=='Figure1':
            # Coordinates must already be projected in kilometres in the same CRS.
            cell=float(d['cell_size_km']);a=np.floor(d['coordinates_a_km']/cell).astype(int);b=np.floor(d['coordinates_b_km']/cell).astype(int)
            cells=np.unique(np.vstack([a,b]),axis=0)
            ca=np.array([(a==c).all(1).sum() for c in cells]);cb=np.array([(b==c).all(1).sum() for c in cells])
            source.update(grid_cells=cells,counts_a=ca,counts_b=cb,shared=(ca>0)&(cb>0),cell_size_km=cell)
            add('overlap','identity_scatter',x=ca,y=cb,xlabel='Cohort A per cell',ylabel='Cohort B per cell')
    elif figure=='Figure3':
        source=relocation_bootstrap(d,repetitions,seed)
        for j in range(2):add('bootstrap'+str(j+1),'bootstrap_histogram',null=source['null_difference'],observed=source['observed'][j],p=source['p'][j],xlabel='Difference in Cohen d')
    elif figure=='ExtendedDataFigure1':
        values=np.atleast_2d(d['measurements']);stats=[anova_icc(values[:,j],d['province']) for j in range(values.shape[1])]
        columns=['n','groups','F','p','ICC1','ICC2'];matrix=np.array([[s[k] for k in columns] for s in stats])
        source=dict(statistics=matrix,columns=np.array(columns),labels=d['measurement_labels'])
        add('reliability','table',matrix=matrix,columns=np.array(columns,dtype=object),labels=d['measurement_labels'])
    elif figure in ['Figure4','ExtendedDataFigure12']:
        effects=np.asarray(d['task_effects']);order=np.asarray(d['best_order'],int)
        if sorted(order.tolist())!=list(range(effects.shape[0])):raise ValueError('best_order must be one zero-based permutation of tasks')
        ordered=effects[order];task_labels=d['task_labels'][order]
        source=dict(task_effects=effects,best_order=order,ordered_task_effects=ordered)
        x=np.arange(len(order));A=np.c_[np.ones(len(order)),x]
        source['slopes']=np.linalg.lstsq(A,ordered,rcond=None)[0][1]
        for roi in np.asarray(d['example_roi'],int).ravel():
            v=ordered[:,roi];rdm=v[:,None]-v[None,:];source['rdm_roi'+str(roi)]=rdm
            heat('RDM_ROI'+str(roi),rdm,task_labels,task_labels)
            add('task_ROI'+str(roi),'task_scatter',x=x,y=v,color=np.array([.6,.45,.65]),labels=task_labels)
        if 'mediation_estimates' in d:
            source.update({k:d[k] for k in ['mediation_estimates','mediation_low','mediation_high','mediation_labels']})
            add('mediation','intervals',estimate=d['mediation_estimates'],low=d['mediation_low'],high=d['mediation_high'],labels=d['mediation_labels'],xlabel='Path coefficient')
    elif figure in ['ExtendedDataFigure2','SupplementaryFigure3']:
        add('trajectory','curves',x=d['components'],y=d['statistics'],labels=d['series_labels'],xlabel='Components removed',ylabel='Statistic')
    elif figure=='ExtendedDataFigure9':
        # Replicates x sampling proportions x series.
        draws=d['correlation_draws'];series=np.stack([np.median(draws,axis=0),np.quantile(draws,.025,axis=0),np.quantile(draws,.975,axis=0)],axis=-1).transpose(1,0,2)
        source.update(series=series)
        add('sensitivity','ribbons',x=d['proportions'],series=series,labels=d['series_labels'])
    elif figure=='ExtendedDataFigure10':
        full=d['full_map'] if 'full_map' in d else d['full_beta']
        fold=d['fold_maps'] if 'fold_maps' in d else d['fold_beta']
        r=np.array([np.corrcoef(full,v)[0,1] for v in fold]);source=dict(correlations=r,unit_labels=d['unit_labels'])
        add('LOO','loo',y=r,labels=d['unit_labels'])
    elif figure=='ExtendedDataFigure11':
        add('ROI_intervals','intervals',estimate=d['estimate'],low=d['low'],high=d['high'],labels=d['roi_labels'],xlabel='Effect')
    elif figure=='ExtendedDataFigure13':
        from map_statistics import spin_test,standardize
        from specialized import heat as correspondence
        source={};D=len(d['dimension_labels']);labels_all=np.r_[d['dimension_labels'],d['gradient_labels']]
        for cohort in ['china','nations']:
            effects=d[cohort+'_dimension_effects'];M=np.c_[effects,d['gradients']];Z=standardize(M)/np.sqrt(len(M)-1);R=Z.T@Z
            ps=np.ones(R.shape);pp=np.full(R.shape,np.nan);sn=np.empty((len(d['spin_indices']),len(R),len(R)));pn=np.full((len(d[cohort+'_dimension_nulls']),len(R),len(R)),np.nan)
            for j in range(len(R)):
                null=standardize(M[:,j][d['spin_indices']].T).T/np.sqrt(len(M)-1)
                sn[:,j,:]=null@Z;ps[j]=(abs(sn[:,j,:])>=abs(R[j])).mean(0)
            for j in range(D):
                null=standardize(d[cohort+'_dimension_nulls'][:,:,j].T).T/np.sqrt(len(M)-1)
                pn[:,j,:]=null@Z;pp[j]=(abs(pn[:,j,:])>=abs(R[j])).mean(0)
            pn[:,D:,:D]=pn[:,:D,D:].transpose(0,2,1)
            pp[:D,:D]=np.maximum(pp[:D,:D],pp[:D,:D].T);pp[D:,:D]=pp[:D,D:].T
            ps=np.maximum(ps,ps.T);combined=np.where(np.isnan(pp),ps,np.maximum(ps,pp));np.fill_diagonal(combined,np.nan)
            source.update({cohort+'_spin_null':sn,cohort+'_permutation_null':pn,cohort+'_maps':M,cohort+'_R':R,cohort+'_p_spin':ps,cohort+'_p_permutation':pp})
            recipe=correspondence(figure,cohort+'_correspondence',R,labels_all,p=combined)
            recipe.update(spin_null=sn,permutation_null=pn,symmetric_spin=True,symmetric_permutation=True,spin_only_cells=True);panels.append(recipe)
        for j,label in enumerate(d['dimension_labels']):
            x=d['china_dimension_effects'][:,j];y=d['nations_dimension_effects'][:,j];r,sn,p=spin_test(x,y,d['spin_indices'])
            add('dimension'+str(j+1),'scatter',x=x,y=y,lh=np.arange(len(x))<len(x)/2,xlabel='China '+str(label),ylabel='Nations '+str(label),stat_null=sn,stat_observed=r,stat_label='r',stats=f'r = {r:.3f}; p = {p:.3f}')
    elif figure=='ExtendedDataFigure15':
        source=joint_pca_axis_angles(d['china_correlations'],d['world_correlations'])
        for k,a in source.items():heat(k,a)
    elif figure=='SupplementaryFigure1':
        x=d['word_scores'];y=d['validity_score'];source=dict(word_scores=x,validity_score=y)
        source['r']=np.array([np.corrcoef(v,y)[0,1] for v in x.T])
        for j in range(x.shape[1]):add('validity'+str(j+1),'scatter',x=x[:,j],y=y,lh=np.ones(len(y),bool),xlabel=str(d['score_labels'][j]),ylabel='Validity score',stats=f'r = {source["r"][j]:.3f}')
    elif figure=='SupplementaryFigure2':
        source=traveling(d['participant1'],d['participant2'])
        source.update(centers=d['center_labels'],common_centers=d['center_labels'],scan_counts=np.array([len(d['center_labels'])]*2))
        for j in range(2):
            heat('deviation'+str(j+1),source['heatmaps'][j],d['center_labels'])
            v=source['relative_error'][j];add('calibration'+str(j+1),'bars',y=v,labels=source['models'],color=np.tile([.5,.6,.75],(len(v),1)),ylabel='Relative error (%)')
    elif figure=='SupplementaryFigure4':
        from collections import Counter
        duration=np.asarray(d['duration']);birth=np.asarray(d['birth_province']).astype(str);dest=np.asarray(d['destination_province']).astype(str)
        if duration.ndim!=1 or not len(duration) or not (len(duration)==len(birth)==len(dest)) or np.any(duration<0) or np.any(duration!=duration.astype(int)):raise ValueError('Aligned nonnegative integer durations and province labels required')
        n=len(duration);routes=list(zip(birth.astype(str),dest.astype(str)));counts=Counter(routes)
        selected=[v for v in sorted(counts,key=lambda k:(-counts[k],k)) if v[0]!=v[1]][:15]
        rc=[]
        for a,b in selected:
            t=duration[(birth==a)&(dest==b)];rc.append([(t<=1).sum(),((t>=2)&(t<=3)).sum(),((t>=4)&(t<=5)).sum(),(t>=6).sum()])
        source=dict(n=n,years=np.arange(int(duration.max())+1),counts=np.bincount(duration.astype(int)),median=np.median(duration),iqr=np.quantile(duration,[.25,.75]),route_labels=np.array([a+' → '+b for a,b in selected]),route_counts=np.asarray(rc).reshape(-1,4))
        source['route_percentages']=100*source['route_counts']/n
        add('duration','bars',y=100*source['counts']/n,labels=source['years'],color=np.tile([.53,.60,.71],(len(source['years']),1)),ylabel='Participants (%)')
        if selected:heat('routes',source['route_percentages'],source['route_labels'],np.array(['0–1','2–3','4–5','6+'],object))
    else:raise ValueError('Unknown figure')
    return source,panels

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--figure',choices=FIGURES,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--repetitions',type=int,default=10000);p.add_argument('--seed',type=int,default=42);a=p.parse_args()
    with np.load(a.input,allow_pickle=False) as f:d=dict(f)
    if a.repetitions<2:raise ValueError('At least two repetitions required')
    source,panels=assemble(a.figure,d,a.repetitions,a.seed)
    print(save(a.figure,a.output,source,panels))
