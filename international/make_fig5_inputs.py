"""Write synthetic individual, unit and spatial files for the Fig. 5 example."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from run_example import synthetic_cohort,NATIONS_BLOCKS,USA_MACROS

def create(folder):
    if folder.exists():raise FileExistsError('Choose a new synthetic input folder')
    nations='Australia|Austria|Belgium|Canada|Colombia|Finland|France|Germany|Italy|Japan|Mexico|Netherlands|Poland|Russia|SouthKorea|Spain|Sweden|Switzerland|United Kingdom|China|Denmark|United States|Greece|Portugal|Hungary|Czechia'.split('|')
    states='Colorado Maryland California Oregon Pennsylvania Missouri Michigan Louisiana Connecticut Wisconsin Georgia New_York Texas New_Jersey Kansas Ohio Virginia Utah'.split()
    rng=np.random.default_rng(22);roi=np.array([('lh_' if i<6 else 'rh_')+f'example{i}' for i in range(12)])
    for cohort,units,seed in [('Nations',nations,101),('USA',states,102)]:
        f=folder/cohort;f.mkdir(parents=True);d=synthetic_cohort(len(units),seed,len(roi));N=len(d['age']);keys=[f'row{j}' for j in range(N)]
        meta=pd.DataFrame(dict(row_key=keys,unit=np.array(units)[d['g']],batch=d['batch'],age=d['age'],female=d['female'],tiv=d['tiv'],age_observed=1,sex_observed=1));meta.to_csv(f/'participants.csv',index=False)
        cortex=pd.DataFrame(d['raw'],columns=roi);cortex.insert(0,'row_key',keys);cortex.iloc[rng.permutation(N)].to_csv(f/'cortex.csv',index=False)
        cols=sum(NATIONS_BLOCKS.values(),[]) if cohort=='Nations' else USA_MACROS;u=pd.DataFrame(d['macro'][:,:len(cols)],columns=cols);u.insert(0,'unit',units)
        if cohort=='USA':u['USA_TL']=45+6*d['culture'][:,0];u['abbreviation']='CO MD CA OR PA MO MI LA CT WI GA NY TX NJ KS OH VA UT'.split()
        else:u['PC1']=d['culture'][:,0];u['PC2']=d['culture'][:,1]
        u.iloc[rng.permutation(len(units))].to_csv(f/'units.csv',index=False);np.savez_compressed(f/'pool.npz',unit_names=np.array(units),pool=d['pool'])
    targets=pd.DataFrame(rng.normal(size=(len(roi),5)),columns=['eigenmode2','SA','China_PC1','China_PC2','China_TL']);targets.insert(0,'roi',roi);targets.iloc[::-1].to_csv(folder/'targets.csv',index=False)
    np.savez_compressed(folder/'spins.npz',roi_names=roi,indices=np.stack([rng.permutation(len(roi)) for _ in range(24)],axis=1))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);create(p.parse_args().output)
