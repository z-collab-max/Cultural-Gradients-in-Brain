"""Rebuild ED5 observed maps and province nulls from individual ComBat inputs."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ.setdefault(key,'1')
from pathlib import Path
import argparse
import h5py
import numpy as np
from scipy.io import savemat

from china_mixed import CONTROLS, read_design, fit_many, province_draws

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--indices',type=Path,default=None)
    parser.add_argument('--protected-pcs',type=int,choices=[2,3],default=2)
    parser.add_argument('--null-protected-pcs',type=int,choices=[2,3])
    parser.add_argument('--draws',type=int,default=10000)
    parser.add_argument('--seed',type=int,default=42)
    args=parser.parse_args()
    if args.null_protected_pcs is None:args.null_protected_pcs=args.protected_pcs
    if args.null_protected_pcs!=args.protected_pcs:raise ValueError('Observed and null ComBat protection must match')
    args.output.mkdir(parents=True,exist_ok=True)
    with h5py.File(args.inputs) as f:
        x,groups=read_design(f['model_design']);birth=f['birth'][()].T
        culture=f['culture'][()].T;yc=f['Y_center'][()].T
        ys={k:f['harmonized'][k][()].T for k in f['harmonized']}
        roi=np.array([''.join(chr(int(c)) for c in f[r][()].ravel()) for r in f['roi'][()].ravel()],dtype=object)
    y=ys[f'PC{args.protected_pcs}_recruit1'];n=len(y)
    fitted=fit_many(x,y,birth[:,:2],groups)
    observed=[(fitted['beta'][pc],fitted['t'][pc],fitted['p'][pc]) for pc in range(2)]
    expcs=fit_many(x,yc,birth[:,:2],groups)['beta'].T
    version='PC12_recruit' if args.protected_pcs==2 else 'PC123_recruit'
    map_name='EDFig5_ComBat_PC12_maps.mat' if args.protected_pcs==2 else 'EDFig5_ComBat_PC123_maps.mat'
    savemat(args.output/map_name,dict(coeff=np.column_stack([v[0] for v in observed]),
        tvalue=np.column_stack([v[1] for v in observed]),pvalue=np.column_stack([v[2] for v in observed]),
        expcs=expcs,brain_var=roi,sample_n=n,protectedPCs=args.protected_pcs,version=version),do_compression=True)
    _,first,group=np.unique(groups[:,0],return_index=True,return_inverse=True)
    indices=None
    if args.indices:
        fixed=np.load(args.indices);indices=fixed['indices_zero_based'][:,:args.draws].T
        if not np.allclose(culture,birth[first,:2]):raise ValueError('Culture rows must match sorted province identifiers')
    predictors=[]
    for pc in range(2):
        px,indices=province_draws(birth[:,pc],groups,args.draws,args.seed,indices)
        predictors.append(px[:,1:])
    models={'Expcs':yc,'CombatMod0':ys['no_culture'],'CombatMod1':y}
    with h5py.File(args.output/'Permutation_Results_FullStats_CulturePC.h5','w') as out:
        out.attrs['protected_birth_pc_count']=args.protected_pcs
        out.attrs['seed']=args.seed
        out.create_dataset('shuffle_indices_zero_based',data=indices)
        for name,y_model in models.items():
            result=fit_many(x,y_model,np.column_stack(predictors),groups)
            dst=out.create_group(f'AllPermRes/{name}')
            for field,key in [('Beta','beta'),('T','t'),('P','p')]:
                dst.create_dataset(field,data=result[key].reshape(2,args.draws,y.shape[1]),compression='gzip')
            print(f'Saved {name}: two PCs, {args.draws} permutations',flush=True)

if __name__=='__main__':main()
