from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analyses'))
import numpy as np
from scipy.io import savemat
from harmonization import harmonize

def prepare(d,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    base=np.c_[d['biological_covariates'],d['culture_pc1'],d['culture_pc2']]
    for name,X in [('protect_PC12',base),('protect_PC12_recruit',np.c_[base,d['recruit']])]:
        y=harmonize(d['cortex'],d['batch'],X)
        savemat(output/(name+'.mat'),dict(cortex=y,protected=X,batch=d['batch']))
        np.savez_compressed(output/(name+'.npz'),cortex=y,protected=X,batch=d['batch'])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    with np.load(a.input,allow_pickle=False) as f:d=dict(f)
    prepare(d,a.output)
