"""Obtain covariate-adjusted parcel measurements from the shared China REML model."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np
from scipy.io import savemat,loadmat

def adjusted_values(raw,covariates,scores,province,city,scanner,pc_index):
    raw=np.asarray(raw,float);covariates=np.asarray(covariates,float);scores=np.asarray(scores,float)
    n=len(raw)
    if raw.ndim!=2 or covariates.shape!=(n,28) or scores.shape!=(n,2):
        raise ValueError('Expected aligned raw parcels, 28 covariates and two cultural scores')
    encoded=[]
    for values in [province,city,scanner]:
        values=np.asarray(values).ravel()
        if len(values)!=n or any(str(v).strip().lower() in ['', 'nan', 'none', '<undefined>'] for v in values):
            raise ValueError('Missing or misaligned province, city or scanner identifiers')
        encoded.append(np.unique(values.astype(str),return_inverse=True)[1])
    if not all(np.isfinite(v).all() for v in [raw,covariates,scores]):raise ValueError('Inputs must be finite')
    pc_index=np.asarray(pc_index,int)
    if pc_index.shape!=(raw.shape[1],) or not np.isin(pc_index,[1,2]).all():raise ValueError('One 1-based cultural PC per parcel is required')
    executable=os.environ.get('MATLAB_EXECUTABLE') or shutil.which('matlab')
    if not executable:raise RuntimeError('MATLAB with Statistics and Machine Learning Toolbox is required; set MATLAB_EXECUTABLE')
    here=Path(__file__).resolve().parent;quote=lambda p:str(p).replace("'","''")
    with tempfile.TemporaryDirectory(prefix='fig2_mixed_') as tmp:
        source=Path(tmp)/'input.mat';target=Path(tmp)/'output.mat'
        savemat(source,dict(raw=raw,covariates=covariates,scores=scores,groups=np.column_stack(encoded[:2]),scanner=encoded[2][:,None],pc_index=pc_index))
        command=f"addpath('{quote(here)}'); fig2_mixed_adjustment('{quote(source)}','{quote(target)}')"
        subprocess.run([executable,'-batch',command],check=True,cwd=here)
        return loadmat(target)['adjusted_mm']
