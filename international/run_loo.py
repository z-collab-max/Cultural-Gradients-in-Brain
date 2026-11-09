"""Export fixed-culture-residual, fold-wise artifact-refitted endpoint maps."""
import argparse
from pathlib import Path
import sys
import numpy as np
from fig5_chain import MODELS,read_cortex,model_inputs
from leave_one_unit_out import leave_one_unit_out
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
from panel_data import save,panel


def run(inputs,output,harmonized,models):
    for model in models:
        cohort=MODELS[model][0]
        if harmonized:
            with np.load(harmonized/(cohort+'_harmonized.npz')) as f: raw=dict(f)
        else:
            raw=read_cortex(inputs,cohort)
        d,*_=model_inputs(raw,inputs,model)
        result=leave_one_unit_out(d)
        result.update(roi_names=raw['roi_names'],unit_names=raw['group_names'])
        if not np.isfinite(result['correlation']).all():raise ValueError('Nonfinite leave-one-out correlations')
        fig='ExtendedDataFigure10'
        save(fig,output/model,result,[panel(fig,model+'_LOO','bars',y=result['correlation'],labels=raw['group_names'].astype(object),ylabel='Spatial correlation')])
        print(model,'LOO',len(result['correlation']),'folds',float(np.min(result['correlation'])),float(np.median(result['correlation'])),float(np.max(result['correlation'])),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['inputs','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--harmonized',type=Path)
    p.add_argument('--models',choices=MODELS,nargs='+',default=list(MODELS))
    a=p.parse_args();run(a.inputs,a.output,a.harmonized,a.models)
