"""Assemble model maps, permutation maps and spatial inputs using a JSON manifest."""
import argparse,json
from pathlib import Path
import numpy as np
from run_specialized import read_input,run

def collect(manifest,output):
    manifest=Path(manifest);spec=json.loads(manifest.read_text(encoding='utf-8'));base=manifest.parent
    d=read_input(base/spec['spatial']);effects=[];nulls=[];names=[];draws=[]
    for model in spec['models']:
        e=read_input(base/model['effect_file']);n=read_input(base/model['null_file'])
        effect=np.asarray(e[model.get('effect_key','beta')]);null=np.asarray(n[model.get('null_key','null_beta')])
        if 'model_index' in model:effect=effect[model['model_index']];null=null[:,model['model_index']]
        effect=effect.reshape(-1)
        if null.ndim!=2 or null.shape[1]!=len(effect):raise ValueError('Null maps must have shape draws x ROI')
        effects.append(effect);nulls.append(null);names.append(model['label'])
        if 'webb_file' in model:
            w=read_input(base/model['webb_file']);v=np.asarray(w[model.get('webb_key','map_draws')])
            if v.ndim!=2 or v.shape[1]!=len(effect):raise ValueError('Webb maps must have shape draws x ROI')
            draws.append(v)
    d.update(effect_maps=np.stack(effects),null_maps=np.stack(nulls,axis=1),model_labels=np.array(names))
    if draws:
        if len(draws)!=len(effects):raise ValueError('Each model needs a Webb file')
        d['webb_map_draws']=np.stack(draws,axis=1)
    return run(spec['analysis'],d,output,spec.get('repetitions',10000),spec.get('seed',42))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();collect(a.manifest,a.output)
