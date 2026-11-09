"""Paired province shuffles of area and volume effect maps in a named ROI order."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
import pandas as pd
from scipy.io import savemat

from china_mixed import CONTROLS, z, read_design, fit_many, province_draws
BASE=CONTROLS

def run(brain,subjects,out,draws,seed):
    out.mkdir(parents=True,exist_ok=True)
    names=json.loads(Path(__file__).with_name('roi_names_307.json').read_text(encoding='utf-8-sig'))
    full_names=json.loads(Path(__file__).with_name('roi_names_308.json').read_text(encoding='utf-8'))
    orders=None
    with h5py.File(brain) as h:
        d=h['brain_data']
        ids=[str(int(x)) for x in d['MRIID'][()].ravel()]
        if len(set(ids))!=len(ids):raise ValueError('Duplicate MRI IDs')
        table=pd.read_csv(subjects,dtype=str).set_index('subject_id',verify_integrity=True)
        provinces=table.loc[ids,'province'].fillna('').str.strip().to_numpy(str)
        def cols(keys):return np.column_stack([d[k][()].ravel() for k in keys])
        for modality in ['area','volume']:
            X,groups=read_design(d)
            selection=np.array([full_names[modality].index(n) for n in names[modality]])
            y=cols(full_names[modality]);c=cols(['pc7_1','pc7_2'])
            center=cols(['centerID_rerange'])[:,0]
            extra=np.empty((len(X),0))
            ok=np.isfinite(np.c_[X,y,c,center,extra]).all(1)
            if np.any(provinces[ok]==''):raise ValueError('Eligible participant with missing province')
            units,first,g=np.unique(groups[ok,0],return_index=True,return_inverse=True)
            if orders is None:
                unit_order=units;rng=np.random.default_rng(seed)
                orders=np.stack([rng.permutation(len(units)) for _ in range(draws)])
            if not np.array_equal(unit_order,units):raise ValueError('Modalities have different province sets')
            scores=c[ok][first]
            if not np.allclose(scores[g],c[ok]):raise ValueError('Cultural scores vary within province')
            _,i=np.unique(center[ok],return_inverse=True);N=np.bincount(i)
            sums=np.zeros((len(N),y.shape[1]));np.add.at(sums,i,y[ok])
            Y=z(y[ok]-(sums/N[:,None])[i])
            for pc in range(2):
                candidates,orders=province_draws(c[ok,pc],groups[ok],draws,seed,orders)
                result=fit_many(X[ok],Y,candidates,groups[ok]);beta=result['beta']
                data=dict(effect=beta[0,selection],t=result['t'][0,selection],p=result['p'][0,selection],df=result['df'][0,selection],spin_source_effect=beta[0],spin_source_roi_names=np.array(full_names[modality]),null_model_effect=beta[0,selection],null=beta[1:,selection],roi_names=np.array(names[modality]),n=len(i),province_names=units,shuffle_indices=orders)
                stem=f'{modality}_PC{pc+1}'
                np.savez_compressed(out/(stem+'.npz'),**data)
                savemat(out/(stem+'.mat'),{k:v.astype(object) if isinstance(v,np.ndarray) and v.dtype.kind=='U' else v for k,v in data.items()},do_compression=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--subjects',type=Path,required=True)
    p.add_argument('--brain',type=Path,required=True);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--draws',type=int,default=10000);p.add_argument('--seed',type=int,default=42)
    a=p.parse_args();run(a.brain,a.subjects,a.output,a.draws,a.seed)

