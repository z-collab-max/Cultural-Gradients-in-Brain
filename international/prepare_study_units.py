"""Assemble named unit tables, comparison targets and spin inputs for the study schema."""
import argparse,json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat
from fig5_chain import MODELS,NATIONS_BLOCKS,USA_MACROS,z,align


def build(config,harmonized,output):
    cfg={k:config.parent/v for k,v in json.loads(config.read_text(encoding='utf-8-sig')).items()}
    def read(key):
        with np.load(cfg[key],allow_pickle=False) as data:return dict(data)
    nations=read('nations_macro');pc1=read('nations_pc1');usa=read('usa_model');um=read('usa_macro')
    joint=loadmat(cfg['nation_culture'],simplify_cells=True)
    scores=pd.DataFrame(joint['score_world'][:,:2],columns=['PC1','PC2']);scores['unit']=joint['matched_nations']
    states=pd.read_excel(cfg['state_culture'],sheet_name='usa_current');states['unit']=states.State.str.replace(' ','_',regex=False)
    roi=np.atleast_1d(np.loadtxt(cfg['roi_order'],dtype=str));targets={'roi':roi}
    for cohort in ['Nations','USA']:
        with np.load(harmonized/(cohort+'_harmonized.npz')) as h:g=h['g'];names=h['group_names'];cortex_roi=h['roi_names']
        if not np.array_equal(cortex_roi,roi):raise ValueError('Reference-map and cortical parcel order differ')
        table=pd.DataFrame(dict(unit=names));sub=output/cohort;sub.mkdir(parents=True,exist_ok=True)
        if cohort=='Nations':
            if not np.array_equal(g,nations['g']) or not np.array_equal(g,pc1['g']):raise ValueError('Country rows differ')
            table[['PC1','PC2']]=align(scores,'unit',names)[['PC1','PC2']].to_numpy()
            for name in sum(NATIONS_BLOCKS.values(),[]):
                vals=pc1['X'][:,list(pc1['labels']).index(name)] if name=='PPPrw' else nations['macro'][:,list(nations['macro_labels']).index(name)]
                grouped=[vals[g==j] for j in range(len(names))]
                if any(np.std(v)>1e-8 for v in grouped):raise ValueError('Nonconstant country macro '+name)
                table[name]=[v.mean() for v in grouped]
            models={'PC1':pc1,'PC2':nations}
            if not np.array_equal(pc1['pool'],nations['pool']):raise ValueError('Separate artifact pools require separate runs')
        else:
            if not np.array_equal(g,usa['g']) or not np.allclose(um['C'],usa['C']):raise ValueError('State rows differ')
            matched=align(states,'unit',names);table['USA_TL']=matched.Cultural_Tightness.to_numpy();table['abbreviation']=matched.Sname.to_numpy()
            for name in USA_MACROS:table[name]=um['macro_group'][:,list(um['labels']).index(name)]
            models={'USA_TL':usa}
        for model,pack in models.items():
            np.testing.assert_allclose(z(table[model].to_numpy()),pack['C'],atol=1e-10,rtol=0)
            for j,key in enumerate(MODELS[model][1:3]):
                values=pack['targets'][:,j]
                if key in targets:np.testing.assert_array_equal(values,targets[key])
                targets[key]=values
        table.to_csv(sub/'units.csv',index=False);pool=next(iter(models.values()))['pool'];np.savez_compressed(sub/'pool.npz',unit_names=names,pool=pool)
    pd.DataFrame(targets).to_csv(output/'targets.csv',index=False);np.savez_compressed(output/'spins.npz',roi_names=roi,indices=np.load(cfg['spin_indices'],allow_pickle=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--harmonized',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();build(a.config,a.harmonized,a.output)
