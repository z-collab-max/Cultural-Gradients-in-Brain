from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'core'))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analyses'))
import numpy as np
from engine import make
from artifact import artifact_endpoint
from protected_reml import fit
from cultural_pca import cultural_pca

def dimension_loo(d,dimensions):
    """Refit the endpoint after omitting each cultural dimension."""
    dimensions=np.asarray(dimensions,float)
    if dimensions.shape[0]!=len(d['C']) or dimensions.shape[1]<3:raise ValueError('Aligned unit-by-dimension matrix with at least three dimensions required')
    a,_=artifact_endpoint(make(d['X'],d['Y'],d),d['C'],d['pool'],'screen_r001');full,_=fit(a,d['C'],a['k']);maps=[];matched=[]
    for omit in range(dimensions.shape[1]):
        scores=cultural_pca(np.delete(dimensions,omit,axis=1))['scores'][:,:2];r=np.array([np.corrcoef(d['C'],v)[0,1] for v in scores.T]);j=int(np.argmax(abs(r)));c=scores[:,j]*np.sign(r[j]);q=dict(d,C=c)
        a,_=artifact_endpoint(make(q['X'],q['Y'],q),c,q['pool'],'screen_r001');beta,_=fit(a,c,a['k']);maps.append(beta);matched.append(j)
    return dict(full_map=full,fold_maps=np.array(maps),matched_components=np.array(matched))
