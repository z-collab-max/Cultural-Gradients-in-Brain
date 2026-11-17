from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'international/core'))
import numpy as np
from combat import combat


def harmonize(cortex,batch,protected):
    # protected contains numeric biological covariates, without an intercept.
    Y=np.asarray(cortex,float);P=np.asarray(protected,float)
    _,b=np.unique(batch,return_inverse=True);B=np.eye(b.max()+1)[b]
    if np.any(B.sum(0)<2):raise ValueError('Each batch needs at least two observations')
    D=np.column_stack([B,P])
    if np.linalg.matrix_rank(D)!=D.shape[1]:raise ValueError('Batch/protected design is rank deficient')
    plan=[dict(ix=np.arange(len(Y)),b=b,B=B,protect=P)]
    return combat(Y,plan)
