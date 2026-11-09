"""Two-intercept REML recovery simulations for SNR and replacement."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np
from scipy.io import savemat,loadmat

def recovery_draws(Y,covariates,scores,province,city,draws):
    Y=np.asarray(Y,float);covariates=np.asarray(covariates,float);scores=np.asarray(scores,float)
    if Y.ndim!=2 or covariates.shape!=(len(Y),28) or scores.shape!=(len(Y),2):
        raise ValueError('Recovery requires aligned cortex, 28 covariates and two cultural scores')
    if draws<3 or not all(np.isfinite(v).all() for v in [Y,covariates,scores]):
        raise ValueError('Use finite inputs and at least three draws per batch')
    groups=[]
    for v in [province,city]:
        v=np.asarray(v).ravel()
        if len(v)!=len(Y) or any(str(k).strip().lower() in ['', 'nan', 'none', '<undefined>'] for k in v):
            raise ValueError('Missing or misaligned birthplace province/city identifiers')
        groups.append(np.unique(v.astype(str),return_inverse=True)[1])
    executable=os.environ.get('MATLAB_EXECUTABLE') or shutil.which('matlab')
    if not executable:raise RuntimeError('MATLAB with Statistics and Machine Learning Toolbox is required; set MATLAB_EXECUTABLE')
    here=Path(__file__).resolve().parent;quote=lambda p:str(p).replace("'","''")
    with tempfile.TemporaryDirectory(prefix='fig2_recovery_') as tmp:
        source=Path(tmp)/'input.mat';target=Path(tmp)/'output.mat'
        savemat(source,dict(Y=Y,covariates=covariates,scores=scores,groups=np.column_stack(groups),draws=draws))
        command=f"addpath('{quote(here)}'); fig2_recovery_mixed('{quote(source)}','{quote(target)}')"
        subprocess.run([executable,'-batch',command],check=True,cwd=here)
        return {k:v for k,v in loadmat(target).items() if not k.startswith('__')}
