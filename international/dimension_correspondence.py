"""Fit six cultural dimensions, China comparisons and gradient correspondence."""
import argparse,json,hashlib,sys
from pathlib import Path
from fig5_chain import read_cortex,model_inputs,make,artifact_endpoint,fit,align,z,unit
from trajectory_fit import fit_levels
import numpy as np
import pandas as pd
from scipy.io import savemat
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflows'));sys.path.insert(0,str(ROOT/'workflows/study_inputs'))
from panel_data import panel,save
from specialized import heat,effect_pca
from matlab_source_reader import read_mat
from trajectory_models import fingerprint


def run(args):
    if args.harmonized:
        with np.load(args.harmonized/'Nations_harmonized.npz') as f:raw=dict(f)
    else:raw=read_cortex(args.inputs,'Nations')
    d,*_=model_inputs(raw,args.inputs,'PC1');roi=d['roi_names'];labels=['PDI','IDV','MAS','UAI','LTO','IVR']
    table=align(pd.read_csv(args.dimensions),'unit',d['group_names']);scores=z(table[labels].to_numpy(float))
    with np.load(args.china_maps) as f:
        china_roi=np.array([s.removesuffix('_thickness') for s in f['roi_names']])
        if set(china_roi)!=set(roi) or len(set(china_roi))!=len(roi):raise ValueError('China atlas labels differ')
        roi_order=np.array([list(china_roi).index(s) for s in roi])
        idx=[list(f['predictor_names']).index(k) for k in ['Hofstede_powerdis_control','Hofstede_indvcoll_control','Hofstede_masvfem_control','Hofstede_unceravoi_control','Hofstede_lonvsho_control','Hofstede_indvres_control']]
        china=f['maps'][idx].T[roi_order]
    with np.load(args.inputs/'spins.npz') as f:
        if not np.array_equal(f['roi_names'],roi):raise ValueError('Rotation ROI order differs')
        spin=f['indices'].T
    refs=np.column_stack([read_mat(args.gradients,k) for k in ['mode2','SAaxis','intervar','xu_evo','reardon_dev','cogPC1']])[roi_order];refs[:,0]*=-1
    names=labels+['eigenmode2','SAAxis','intersubjvar','evoexp','devexp','cogpc1']
    maps=[];nulls=[];orders=np.array([np.random.RandomState(args.seed+i).permutation(len(d['C'])) for i in range(args.draws)])
    for j,label in enumerate(labels):
        if args.dimension and label!=args.dimension:continue
        q=dict(d,C=scores[:,j]);s=make(q['X'],q['Y'],q);a,_=artifact_endpoint(s,q['C'],q['pool'],'screen_r001');b=fit_levels(a,q['C'],[a['k']])[0]
        np.testing.assert_allclose(b,fit(a,q['C'],a['k'])[0],atol=1e-7,rtol=1e-6)
        folder=args.output/'cache'/label;folder.mkdir(parents=True,exist_ok=True)
        signature=dict(inputs=fingerprint(q),draws=args.draws,seed=args.seed,algorithm='endpoint_unrestricted_v1')
        cfg=folder/'configuration.json'
        if cfg.exists() and json.loads(cfg.read_text())!=signature:raise ValueError('Output cache has different inputs')
        cfg.write_text(json.dumps(signature,indent=2));blocks=[]
        for first in range(0,args.draws,100):
            f=folder/f'block_{first:05d}.npz'
            if f.exists():
                with np.load(f) as old:bb=old['beta']
            else:
                if args.assemble_only:raise FileNotFoundError(f)
                bb=[]
                for order in orders[first:first+100]:
                    c=q['C'][order];ap,_=artifact_endpoint(s,c,q['pool'],'screen_r001');bb.append(fit_levels(ap,c,[ap['k']])[0])
                bb=np.array(bb);temporary=f.with_suffix('.tmp.npz');np.savez_compressed(temporary,beta=bb);temporary.replace(f)
            blocks.append(bb)
            if first%1000==0:print(label,first+len(bb),'/',args.draws,flush=True)
        nb=np.concatenate(blocks);maps.append(b);nulls.append(nb);target=china[:,j]
        r=float(np.corrcoef(target,b)[0,1]);pn=unit(nb.T).T@unit(target[:,None]);pn=pn.ravel();sn=unit(target[spin.T]).T@unit(b[:,None]);sn=sn.ravel()
        pp=(abs(pn)>=abs(r)).mean();ps=(abs(sn)>=abs(r)).mean()
        data=dict(roi_names=roi,target=target,culture=b,r=r,permutation_null=pn,spin_null=sn,p_permutation=pp,p_spin=ps,endpoint_k=a['k'])
        panels=[panel('ExtendedDataFigure13',label+'_China','scatter',x=z(target),y=z(b),lh=np.char.startswith(roi,'lh_'),xlabel='China '+label+' effect',ylabel='Nations '+label+' effect',stat_null=sn,stat_observed=r,stat_label='r',stats=f'r = {r:.3f}; p = {ps:.3f}'),
                panel('ExtendedDataFigure13',label+'_shuffle','histogram',null=pn,observed=r,p=pp,xlabel='Spatial correlation')]
        dest=save('ExtendedDataFigure13',args.output/label,data,panels)
        savemat(dest/'surface.mat',dict(roi_names=roi.astype(object),values=z(b),color_limits=[-4,4]))
    if args.dimension:return
    maps=np.array(maps).T;M=np.c_[maps,refs];Z=unit(M);corr=Z.T@Z;ps=np.ones((12,12));pp=np.full((12,12),np.nan)
    spin_null=np.empty((spin.shape[0],12,12));permutation_null=np.full((args.draws,12,12),np.nan)
    for j in range(12):
        sn=unit(M[:,j][spin.T]).T@Z;spin_null[:,j,:]=sn;ps[j]=(abs(sn)>=abs(corr[j])).mean(0)
    for j in range(6):
        pn=unit(nulls[j].T).T@Z;permutation_null[:,j,:]=pn;pp[j]=(abs(pn)>=abs(corr[j])).mean(0)
    for i in range(6,12):pp[i,:6]=pp[:6,i];permutation_null[:,i,:6]=permutation_null[:,:6,i]
    pp[:6,:6]=np.maximum(pp[:6,:6],pp[:6,:6].T)
    ps=np.maximum(ps,ps.T)
    combined=np.where(np.isnan(pp),ps,np.maximum(ps,pp));np.fill_diagonal(combined,np.nan)
    recipe=heat('ExtendedDataFigure13','Nations_dimensions',corr,np.array(names,object),p=combined)
    recipe.update(spin_null=spin_null,permutation_null=permutation_null,symmetric_spin=True,symmetric_permutation=True,spin_only_cells=True)
    save('ExtendedDataFigure13',args.output/'correspondence',dict(spin_null=spin_null,permutation_null=permutation_null,R=corr,p_spin=ps,p_permutation=pp,p_combined=combined,labels=np.array(names),maps=M,roi_names=roi),[recipe])
    culture=[]
    for model in ['PC1','PC2']:
        q,*_=model_inputs(raw,args.inputs,model);a,_=artifact_endpoint(make(q['X'],q['Y'],q),q['C'],q['pool'],'screen_r001');culture.append(fit_levels(a,q['C'],[a['k']])[0])
    data,panels,_=effect_pca(dict(dimension_effects=maps,culture_effects=np.array(culture).T,spin_indices=spin,reference_component=[0,0],dimension_labels=np.array(labels,object),left_hemisphere=np.char.startswith(roi,'lh_')))
    data.update(dimension_effects=maps,culture_effects=np.array(culture).T,roi_names=roi,dimension_labels=np.array(labels))
    panels.append(panel('ExtendedDataFigure15','Nations_effect_biplot','pca_biplot',scores=data['scores'][:,:2],loadings=data['loadings'][:,:2],explained=data['explained'][:2]*100,labels=np.array(labels,object)))
    panels.append(panel('ExtendedDataFigure15','Nations_effect_scree','bars',y=data['explained']*100,labels=np.array([f'PC{j+1}' for j in range(6)],object),ylabel='Variance (%)'))
    dest=save('ExtendedDataFigure15',args.output,data,panels)
    for j in range(2):
        savemat(dest/f'Effect_PC{j+1}_surface.mat',dict(roi_names=roi.astype(object),values=z(data['scores'][:,j]),color_limits=[-4,4]))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['inputs','dimensions','china-maps','gradients','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--harmonized',type=Path);p.add_argument('--draws',type=int,default=10000);p.add_argument('--seed',type=int,default=42)
    p.add_argument('--dimension',choices=['PDI','IDV','MAS','UAI','LTO','IVR']);p.add_argument('--assemble-only',action='store_true')
    args=p.parse_args()
    if args.draws<2:p.error('At least two shuffles are required')
    run(args)
