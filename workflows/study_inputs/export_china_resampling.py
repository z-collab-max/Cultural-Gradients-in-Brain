"""Export China resampling summaries and draws from a generated result file."""
from pathlib import Path
import argparse,sys,numpy as np
from scipy.io import savemat
from matlab_source_reader import read_mat

def generate(result_file,output_dir):
    result_file=Path(result_file);out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    ratios=np.atleast_1d(read_mat(result_file,'ratios'));summary={'ratios':ratios}
    draws={'ratios':ratios,'n_draws':int(read_mat(result_file,'n_perms'))}
    for pc in [1,2]:
        draws[f'map_target_PC{pc}']=read_mat(result_file,f'map_target_{pc}')
        for mode in ['Prov','Subj']:
            key=f'{mode}_PC{pc}';a=np.asarray(read_mat(result_file,f'Res_{mode}_Corr{pc}')).reshape(len(ratios),draws['n_draws'])
            summary[key+'_valid_draws']=np.isfinite(a).sum(axis=1)
            summary[key]=np.c_[np.nanmean(a,axis=1),np.nanquantile(a,[.025,.975],axis=1,method='hazen').T]
            draws[key+'_correlation_draws']=a
            for i,b in enumerate(read_mat(result_file,f'Beta_{mode}_PC{pc}')):draws[f'{key}_beta_ratio{i+1}']=b
    for name,data in [('EDFig9_China',summary),('EDFig9_China_draws',draws)]:
        savemat(out/(name+'.mat'),data,do_compression=True);np.savez_compressed(out/(name+'.npz'),**data)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--result-file',required=True);p.add_argument('--output-dir',required=True)
    a=p.parse_args();generate(a.result_file,a.output_dir)
