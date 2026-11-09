"""Fit variance-threshold maps using raw cortex and the Figure 5 unit design."""
from pathlib import Path
import argparse
import sys
from fig5_chain import MODELS, read_cortex, model_inputs, make, artifact_endpoint, fit
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
from panel_data import panel, save


def threshold_maps(d, model, ref='endpoint'):
    a,_=artifact_endpoint(make(d['X'],d['Y'],d),d['C'],d['pool'],'screen_r001')
    thresholds=np.array([80,85,90,95])
    ks=np.searchsorted(a['cum_explained'],thresholds)+1
    reference_k=a['k'] if ref=='endpoint' else np.searchsorted(a['cum_explained'],99.9)+1
    by_k={k:fit(a,d['C'],int(k))[0] for k in np.unique(np.r_[ks,reference_k])}
    base=by_k[reference_k]; maps=np.column_stack([by_k[k] for k in ks])
    r=np.array([np.corrcoef(base,col)[0,1] for col in maps.T])
    data=dict(x=base,y=maps,thresholds=thresholds,ks=ks,reference_k=reference_k,
              reference_achieved=a['cum_explained'][reference_k-1],achieved=a['cum_explained'][ks-1],correlations=r)
    p=panel('ExtendedDataFigure2',model+'_thresholds','curves',x=thresholds,y=r[:,None],
            labels=np.array([ref],object),xlabel='Explained variance (%)',ylabel='Correlation with reference')
    return data,[p]


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--models',choices=MODELS,nargs='+',default=list(MODELS))
    parser.add_argument('--harmonized',type=Path,help='Optional Figure 5 harmonized arrays directory')
    parser.add_argument('--reference',choices=['endpoint','99.9'],default='endpoint')
    args=parser.parse_args();raw={}
    for model in args.models:
        cohort=MODELS[model][0]
        if cohort not in raw:
            if args.harmonized:
                with np.load(args.harmonized/(cohort+'_harmonized.npz')) as f:raw[cohort]=dict(f)
            else:
                raw[cohort]=read_cortex(args.inputs,cohort)
        d,*_=model_inputs(raw[cohort],args.inputs,model)
        source,panels=threshold_maps(d,model,args.reference)
        save('ExtendedDataFigure2',args.output/model,source,panels)

