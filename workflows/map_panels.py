import numpy as np
from panel_data import panel,labels
from map_statistics import standardize
from webb_statistics import cluster_score_intervals,spatial_statistics


def assemble_maps(figure,data,repetitions,seed,with_ci=False):
    effects=np.atleast_2d(data['effect_maps']);G=standardize(data['gradients']);null=data['null_maps']
    rawspins=np.asarray(data['spin_indices']);spins=rawspins.astype(int);R=G.shape[0];K=G.shape[1]
    if not np.array_equal(rawspins,spins):raise ValueError('Spin indices must be integers')
    if effects.shape[1]!=R or null.shape[1:]!=effects.shape:raise ValueError('Map dimensions disagree')
    if spins.ndim!=2 or spins.shape[1]!=R or np.any(spins<0) or np.any(spins>=R):raise ValueError('Invalid spin indices')
    projector=np.linalg.solve(G.T@G,G.T);names=np.asarray(data.get('gradient_labels',labels('Gradient',K)),object)
    panels=[];source=dict(effect_maps=effects,gradients=G,gradient_labels=names)
    for m,effect in enumerate(effects):
        name=f'{figure}_Model{m+1}';point=standardize(effect);beta=projector@point
        B=standardize(null[:,m,:].T).T;nb=B@projector.T
        rotated=point[spins];rotated=standardize(rotated.T).T;sb=rotated@projector.T
        pperm=(abs(nb)>=abs(beta)).mean(0)
        pspin=(abs(sb)>=abs(beta)).mean(0)
        source.update({f'model{m+1}_beta':beta,f'model{m+1}_null_beta':nb,f'model{m+1}_spin_beta':sb,
                       f'model{m+1}_p_permutation':pperm,f'model{m+1}_p_spin':pspin})
        panels.append(panel(figure,name+'_gradients','violin',beta=beta,null=nb,additional_null=sb,additional_observed=beta,p=np.maximum(pperm,pspin),labels=names))
        for k in range(min(2,K)):
            panels.append(panel(figure,name+f'_scatter{k+1}','scatter',x=G[:,k],y=point,
                                lh=np.asarray(data.get('left_hemisphere',np.arange(R)<R//2),bool),
                                stat_null=sb[:,k],stat_observed=beta[k],stats=fr'\beta = {beta[k]:.3f}, p_{{spin}} = {pspin[k]:.3f}',xlabel=str(names[k]),ylabel='Cultural effect'))
            panels.append(panel(figure,name+f'_null{k+1}','histogram',null=nb[:,k],observed=beta[k],p=pperm[k],xlabel='Partial coefficient'))
        if with_ci:
            if figure=='Figure5':
                if 'webb_map_draws' not in data:raise ValueError('Figure5 requires Webb map draws from international endpoint refits')
                ci=spatial_statistics(effect,data['webb_map_draws'][:,m,:],G,quantile_method='linear')
            elif 'cluster_influence' in data:
                ci=cluster_score_intervals(effect,data['cluster_influence'][m],G,repetitions,seed+m)
            elif 'webb_map_draws' in data:
                ci=spatial_statistics(effect,data['webb_map_draws'][:,m,:],G,quantile_method='linear')
            else:raise ValueError('CI requires cluster_influence or Webb map draws, not permutation nulls')
            source.update({f'model{m+1}_webb_{k}':v for k,v in ci.items()})
            panels.append(panel(figure,name+'_Webb','intervals',estimate=ci['partial_beta'],low=ci['beta_ci'][0],high=ci['beta_ci'][1],labels=names,xlabel='Partial coefficient'))
    return source,panels
