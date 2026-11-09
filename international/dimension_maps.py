"""Refit international cultural dimensions and leave-one-dimension-out PC maps."""
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from fig5_chain import align, model_inputs, read_cortex, make, artifact_endpoint, fit, z
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
from specialized import effect_pca, heat
from panel_data import save
from cultural_pca import cultural_pca
from map_statistics import spin_test


def fitted(d,c):
    a,_=artifact_endpoint(make(d['X'],d['Y'],dict(d,C=c)),c,d['pool'],'screen_r001')
    return fit(a,c,a['k'])[0]


def run(inputs,dimensions,out,components,harmonized):
    if harmonized:
        with np.load(harmonized/'Nations_harmonized.npz') as f:raw=dict(f)
    else:
        raw=read_cortex(inputs,'Nations')
    if dimensions.suffix.lower()=='.xlsx':
        world=pd.read_excel(dimensions,sheet_name='Culdis_world')
        world['unit']=world.region.replace({'Czech Republic':'Czechia'})
        world=world[world.unit.isin(raw['group_names'])]
        fields=['PDI','IDV','MAS','UAI','LTO','IVR']
        table=align(world,'unit',raw['group_names'])[[s+'_cc' for s in fields]]
        table.columns=fields
    else:
        table=align(pd.read_csv(dimensions),'unit',raw['group_names'])
    with np.load(inputs/'spins.npz') as f:
        if not np.array_equal(f['roi_names'],raw['roi_names']):raise ValueError('Spin ROI order differs from cortical data')
        spin=f['indices'].T
    dims=table.to_numpy(float);labels=np.asarray(table.columns,str)
    models=[model_inputs(raw,inputs,m)[0] for m in ['PC1','PC2']]
    original=np.stack([d['C'] for d in models])
    full=np.stack([fitted(d,d['C']) for d in models])
    maps=np.column_stack([fitted(models[0],z(col)) for col in dims.T])
    folds=[];matching=[];cor=[];pv=[]
    for omit in range(len(labels)):
        scores=cultural_pca(np.delete(dims,omit,axis=1))['scores'][:,:2]
        r=np.corrcoef(original,scores.T)[:2,2:]
        order=np.array([0,1]) if abs(np.diag(r)).sum()>=abs(np.diag(r[:,::-1])).sum() else np.array([1,0])
        row=[fitted(d,z(scores[:,k])*np.sign(r[j,k])) for j,(d,k) in enumerate(zip(models,order))]
        tests=[spin_test(full[j],v,spin) for j,v in enumerate(row)]
        folds.append(row);matching.append(order);cor.append([t[0] for t in tests]);pv.append([t[2] for t in tests])
    source=dict(full_maps=full,fold_maps=np.array(folds),matching=np.array(matching),correlations=np.array(cor),p_spin=np.array(pv),dimension_labels=labels,roi_names=raw['roi_names'])
    save('CulturalDimensionOmission',out,source,[heat('CulturalDimensionOmission','Dimension_omission',source['correlations'],labels,['PC1','PC2'],source['p_spin'])])
    d=dict(dimension_effects=maps,culture_effects=full.T,spin_indices=spin,reference_component=components,dimension_labels=labels,left_hemisphere=np.char.startswith(raw['roi_names'],'lh_'))
    data,panels,_=effect_pca(d);data.update(dimension_effects=maps,culture_effects=full.T,roi_names=raw['roi_names'],dimension_labels=labels)
    save('ExtendedDataFigure15',out,data,panels)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['inputs','dimensions','output']:p.add_argument('--'+key,required=True,type=Path)
    p.add_argument('--components',required=True,type=int,nargs=2,help='Two zero-based effect PC indices')
    p.add_argument('--harmonized',type=Path)
    a=p.parse_args();run(a.inputs,a.dimensions,a.output,a.components,a.harmonized)

