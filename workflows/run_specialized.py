"""Build figure-specific source files and a MATLAB rendering manifest."""
from pathlib import Path
import argparse,json,re,sys
import numpy as np
from scipy.io import loadmat,savemat
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analyses'))
from specialized import *
from panel_data import save
from map_panels import assemble_maps

FUNCTIONS={'province_summary':province_summary,'demographics':demographics,'gradients':gradients,
 'commonality':commonality_panels,'hierarchy':hierarchy,'rice_culture':correlation_table,
 'effect_pca':effect_pca,'resampling':resampling,'word_layout':word_layout,'traveling_unshared':traveling_unshared,
 'ecological_comparison':ecological_comparison,'psychometrics':psychometric_table,'brain_icc':brain_icc}
FIGURES={'province_summary':'Figure2','demographics':'ExtendedDataFigure1','gradients':'ExtendedDataFigure3',
 'commonality':'Figure3','hierarchy':'Figure4','rice_culture':'ExtendedDataFigure8','effect_pca':'ExtendedDataFigure15',
 'resampling':'ExtendedDataFigure9','word_layout':'SupplementaryFigure1','traveling_unshared':'SupplementaryFigure2',
 'network':'ExtendedDataFigure11','combat_versions':'ExtendedDataFigure5','area_volume':'ExtendedDataFigure7',
 'rice_brain':'ExtendedDataFigure8','international':'Figure5','dimension_omission':'CulturalDimensionOmission',
 'reference_histograms':'ExtendedDataFigure6','thresholds':'ExtendedDataFigure2','overlap_curves':'ExtendedDataFigure9',
 'ecological_comparison':'ExtendedDataFigure1','psychometrics':'ExtendedDataFigure1','brain_icc':'ExtendedDataFigure1'}

def read_input(path):
    path=Path(path)
    if path.suffix=='.mat':return {k:v for k,v in loadmat(path,simplify_cells=True).items() if not k.startswith('__')}
    with np.load(path,allow_pickle=False) as f:return dict(f)

def run(analysis,d,output,repetitions=10000,seed=42):
    if analysis not in FIGURES:raise ValueError('Unknown analysis '+analysis)
    if repetitions<2:raise ValueError('At least two resampling draws required')
    for k,v in d.items():
        a=np.asarray(v)
        if a.dtype.kind in 'fc' and not np.isfinite(a).all():raise ValueError('Nonfinite input '+k)
    figure=FIGURES[analysis];layout=None
    if analysis in FUNCTIONS:source,panels,layout=FUNCTIONS[analysis](d)
    elif analysis=='network':source,panels,layout=network(d,repetitions,seed)
    elif analysis in ['combat_versions','area_volume','rice_brain','international']:
        source,panels=assemble_maps(figure,d,repetitions,seed,with_ci=analysis=='international')
        if 'model_labels' not in d or len(d['model_labels'])!=len(np.atleast_2d(d['effect_maps'])):raise ValueError('One explicit model label per effect map is required')
        source['model_labels']=np.asarray(d['model_labels']);source['null_maps']=d['null_maps'];source['spin_indices']=d['spin_indices']
        for j in range(len(d['model_labels'])):
            for p in panels:
                if f'Model{j+1}_' in p['name'] and 'ylabel' in p:p['ylabel']=str(d['model_labels'][j])
    elif analysis=='reference_histograms':
        maps=np.atleast_2d(d['effect_maps']);null=np.asarray(d['null_maps'])
        refs=np.broadcast_to(d['reference'],maps.shape)
        obs=np.array([pearsonr(v,ref).statistic for v,ref in zip(maps,refs)])
        null_r=np.empty(null.shape[:2])
        for j,ref in enumerate(refs):
            null_r[:,j]=standardize(null[:,j,:].T).T @ standardize(ref)/(len(ref)-1)
        p=(abs(null_r)>=abs(obs)).mean(0);source=dict(observed=obs,null_r=null_r,p=p);panels=[]
        for j,r in enumerate(obs):panels.append(panel(figure,f'Reference_histogram{j+1}','histogram',null=null_r[:,j],observed=r,p=p[j],xlabel='Spatial correlation'))
    elif analysis=='dimension_omission':
        full=np.asarray(d['full_maps']);fold=np.asarray(d['fold_maps']);spins=d['spin_indices'];cor=np.empty((len(fold),len(full)));pv=cor.copy();aligned=fold.copy();matching=[];cell_nulls=np.empty((len(spins),len(fold),len(full)))
        for i,f in enumerate(fold):
            # Fold maps must already be aligned using cultural scores upstream.
            order=tuple(range(len(full)))
            matching.append(order)
            for j,k in enumerate(order):
                aligned[i,j]=f[k];cor[i,j],cell_nulls[:,i,j],pv[i,j]=spin_test(full[j],aligned[i,j],spins)
        source=dict(spin_null=cell_nulls,correlations=cor,p_spin=pv,aligned_maps=aligned,matching=np.array(matching),dimension_labels=d['dimension_labels']);panels=[dict(heat(figure,'Dimension_omission_correlations',cor,d['dimension_labels'],labels('PC',len(full)),pv),cell_nulls=cell_nulls)]
    elif analysis=='thresholds':
        maps=np.asarray(d['threshold_maps']);base=np.asarray(d['endpoint_map']);r=np.array([np.corrcoef(base,x)[0,1] for x in maps]);source=dict(correlations=r,threshold_labels=d['threshold_labels']);panels=[panel(figure,'Endpoint_threshold_comparison','curves',x=np.arange(len(r)),y=r[:,None],labels=np.array(['Endpoint'],object),xlabel='Threshold index',ylabel='Correlation with endpoint')]
    elif analysis=='overlap_curves':
        a=np.asarray(d['overlap_counts']);total=a.sum(1)
        if np.any(total<=0):raise ValueError('Each scale must contain occupied cells')
        source=dict(scales=d['scales'],counts=a,percentages=100*a/total[:,None]);panels=[panel(figure,'Grid_overlap_curve','curves',x=d['scales'],y=source['percentages'],labels=np.array(['Both','Imaging only','Culture only'],object),xlabel='Grid size (km)',ylabel='Occupied cells (%)')]
    folder=Path(output)/analysis;folder.mkdir(parents=True,exist_ok=True)
    for p in panels:p['name']=analysis+'_'+p['name']
    dest=save(figure,folder,source,panels);jobs=[]
    if analysis=='word_layout':
        matlab_source=dict(source)
        for key in ['word_labels','dictionary']:
            cells=np.empty(np.asarray(source[key]).shape,object)
            for index in np.ndindex(cells.shape):cells[index]=np.frombuffer(str(source[key][index]).encode('utf-16-le'),dtype='<u2')
            matlab_source[key]=cells
        savemat(dest/'source_data.mat',matlab_source,do_compression=True)
    if panels:jobs.append(dict(kind='native',input='panels.mat'))
    if layout:jobs.append(dict(kind='layout',function=layout,input='source_data.mat'))
    if 'roi_names' in d:
        maps=[]
        if analysis=='gradients':maps=[(f'gradient{j+1}',source['maps'][:,j]) for j in range(source['maps'].shape[1])]
        elif analysis=='hierarchy':maps=[(f'branch{j+1}_{metric}',source[f'branch{j+1}_{metric}']) for j in range(len(d['branch_labels'])) for metric in ['rsa_map','slope']]
        elif analysis=='effect_pca':maps=[(f'effect_pc{j+1}',source['scores'][:,j]) for j in range(2)]
        elif 'effect_maps' in d:maps=[(f'model{j+1}',v) for j,v in enumerate(np.atleast_2d(d['effect_maps']))]
        for name,v in maps:
            limit=max(float(np.max(abs(v))),1e-12);filename=name+'_surface.mat';savemat(dest/filename,dict(roi_names=np.asarray(d['roi_names'],object),values=v,color_limits=[-limit,limit]));jobs.append(dict(kind='surface',input=filename,stem=name))
    (dest/'render_manifest.json').write_text(json.dumps(jobs,indent=2),encoding='utf-8')
    print(dest);return dest

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--analysis',choices=sorted(FIGURES),required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--repetitions',type=int,default=10000);p.add_argument('--seed',type=int,default=42);a=p.parse_args()
    run(a.analysis,read_input(a.input),a.output,a.repetitions,a.seed)
if __name__=='__main__':main()
