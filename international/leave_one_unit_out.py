from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'core'))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analyses'))
import numpy as np
from engine import make
from core import state,unit
from artifact import artifact_endpoint
from protected_reml import fit,_legacy_fit
from loo_protection import directions,clean


def leave_one_unit_out(d):
    full,_=artifact_endpoint(make(d['X'],d['Y'],d),d['C'],d['pool'],'screen_r001')
    full_beta,_,_=fit(full,d['C'],full['k'],True)
    _,anchor=_legacy_fit(full,d['C'],0);E,_=directions(full,d['C'],anchor)
    maps=[];correlations=[];counts=[]
    for omit in range(len(d['C'])):
        keep=d['g']!=omit;groups=np.delete(np.arange(len(d['C'])),omit)
        q=dict(d,g=np.searchsorted(groups,d['g'][keep]),C=d['C'][groups],X=d['X'][keep],Y=d['Y'][keep])
        q['pool']=np.array([np.searchsorted(groups,p[p!=omit]) for p in d['pool']])
        a,_=artifact_endpoint(make(q['X'],q['Y'],q),q['C'],q['pool'],'screen_r001')
        _,lam=_legacy_fit(a,q['C'],0)
        e=E[keep];w=lam[None,:]/(1+a['n'][:,None]*lam[None,:])
        H=e-w[a['g']]*(a['Z'].T@e)[a['g']]
        y,_,_=clean(a,H)
        beta,_=_legacy_fit(state(q['X'],y,q['g'],q['C'],q['targets']),q['C'],0)
        maps.append(beta);correlations.append(np.corrcoef(beta,full_beta)[0,1]);counts.append(a['k'])
    return dict(full_beta=full_beta,fold_beta=np.array(maps),correlation=np.array(correlations),artifact_pc_count=np.array(counts))
