import numpy as np
from map_statistics import standardize


def cultural_pca(values):
    X=standardize(values);u,s,v=np.linalg.svd(X,full_matrices=False)
    signs=np.sign(v[np.arange(len(v)),np.argmax(abs(v),axis=1)])
    return dict(scores=u*s*signs,loadings=v.T*signs,explained=s*s/(s@s),
                mean=np.mean(values,axis=0),scale=np.std(values,axis=0,ddof=1))


def leave_one_dimension_out(values):
    full=cultural_pca(values);rows=[]
    for j in range(values.shape[1]):
        fit=cultural_pca(np.delete(values,j,axis=1));a=full['scores'][:,:2];b=fit['scores'][:,:2]
        corr=np.corrcoef(a.T,b.T)[:2,2:]
        rows.append(dict(omitted=j,scores=fit['scores'],cross_correlations=corr))
    return full,rows


def axis_angles(loadings_a,loadings_b):
    a=np.asarray(loadings_a);b=np.asarray(loadings_b)
    cosine=(a.T@b)/(np.linalg.norm(a,axis=0)[:,None]*np.linalg.norm(b,axis=0)[None,:])
    return np.degrees(np.arccos(np.clip(abs(cosine),0,1)))


def joint_pca_axis_angles(correlations_china,correlations_world):
    # Rows are joint PCs; columns are local PCs. Use the displayed two-PC plane.
    a=np.asarray(correlations_china)[:2,:2];b=np.asarray(correlations_world)[:2,:2]
    return dict(china_vs_joint=axis_angles(a,np.eye(2)),world_vs_joint=axis_angles(b,np.eye(2)),china_vs_world=axis_angles(a,b))


def cronbach_alpha(items):
    items=np.asarray(items,float);k=items.shape[1]
    if k<2 or not np.isfinite(items).all():raise ValueError('Complete item matrix required')
    return k/(k-1)*(1-items.var(0,ddof=1).sum()/items.sum(1).var(ddof=1))


def joint_cultural_pca(provinces,nations,china_sample_mean,china_world_score,anchor_dimension=5):
    """Align province scores to the China world reference, then fit three PCAs."""
    provinces=np.asarray(provinces,float);nations=np.asarray(nations,float)
    offset=np.asarray(china_world_score)-np.asarray(china_sample_mean)
    if provinces.shape[1]!=nations.shape[1] or offset.shape!=(provinces.shape[1],):raise ValueError('Cultural dimensions must align')
    joint=cultural_pca(np.vstack([provinces+offset,nations]));local=cultural_pca(provinces);world=cultural_pca(nations)
    for model in [joint,local]:
        for j in range(2):
            if model['loadings'][anchor_dimension,j]<0:model['loadings'][:,j]*=-1;model['scores'][:,j]*=-1
    n=len(provinces)
    for j in range(2):
        if np.corrcoef(joint['scores'][n:,j],world['scores'][:,j])[0,1]<0:world['loadings'][:,j]*=-1;world['scores'][:,j]*=-1
    k=min(3,joint['scores'].shape[1]);a=joint['scores'][:n,:k];b=joint['scores'][n:,:k]
    return dict(score_joint_provinces=joint['scores'][:n],score_joint_nations=joint['scores'][n:],
                coeff_joint=joint['loadings'],explained_joint=100*joint['explained'],
                corr_matrix=np.corrcoef(a.T,local['scores'][:,:k].T)[:k,k:],
                corr_matrix_world=np.corrcoef(b.T,world['scores'][:,:k].T)[:k,k:],
                local_scores=local['scores'],world_scores=world['scores'],alignment_offset=offset)
