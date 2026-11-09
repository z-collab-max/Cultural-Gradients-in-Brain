import numpy as np


def remove_artifacts(cortex,score_pool,coefficient_pool,pool_correlations,n_components):
    # score_pool: individuals x shuffles; coefficient_pool: shuffles x ROIs.
    Y=np.asarray(cortex,float);S=np.asarray(score_pool,float);B=np.asarray(coefficient_pool,float)
    keep=abs(np.asarray(pool_correlations).ravel())<.01
    if S.shape!=(len(Y),len(keep)) or B.shape!=(len(keep),Y.shape[1]):
        raise ValueError('Incompatible cortex, score and coefficient dimensions')
    if not 0<=n_components<min(len(Y),keep.sum()):raise ValueError('Invalid component count')
    cleaned=np.empty_like(Y);variance=[]
    for j in range(Y.shape[1]):
        contribution=S[:,keep]*B[keep,j]
        contribution-=contribution.mean(0)
        u,s,_=np.linalg.svd(contribution,full_matrices=False)
        design=np.column_stack([np.ones(len(Y)),u[:,:n_components]*s[:n_components]])
        cleaned[:,j]=Y[:,j]-design@np.linalg.lstsq(design,Y[:,j],rcond=None)[0]
        variance.append(s*s/(s@s))
    return cleaned,np.array(variance)
