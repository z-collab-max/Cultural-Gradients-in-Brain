from pathlib import Path
import numpy as np
import pandas as pd
from cultural_pca import cultural_pca,leave_one_dimension_out,axis_angles,cronbach_alpha
from psychometrics import anova_icc,rwg
from group_summaries import summaries,shared_roi_bootstrap,grid_overlap
from map_statistics import map_correlation,commonality,partial_map_coefficients,bh_fdr
from traveling_calibration import fit_predict
from matrix_io import save_source


def run(output):
    rng=np.random.default_rng(20261103);groups=np.repeat(np.arange(30),20)
    items=np.clip(np.rint(3+rng.normal(size=(600,6))),1,5)
    psych=anova_icc(items.mean(1),groups);assert psych['n']==600
    frame=pd.DataFrame(items,columns=list('abcdef'));frame['from_province']=groups
    agreement,_=rwg(frame,list('abcdef'),5);alpha=cronbach_alpha(items)
    country=rng.normal(size=(30,7));pca,loo=leave_one_dimension_out(country)
    assert len(loo)==7 and np.isclose(pca['explained'].sum(),1)
    assert np.max(abs(np.diag(axis_angles(pca['loadings'],pca['loadings']))))<1e-5
    save_source(output,'cultural_pca',**pca)
    save_source(output,'psychometrics',alpha=alpha,rwg=agreement,**psych)
    summary=summaries(items,groups);save_source(output,'province_summary',**summary)
    boot=shared_roi_bootstrap(items[:,0],items[:,1],groups,20);save_source(output,'roi_bootstrap',**boot)
    components=rng.normal(size=(100,3));target=components@np.array([.3,.4,.2])+rng.normal(size=100)
    parts,coef,pred=commonality(target,*components.T);assert np.isclose(parts[:4].sum(),parts[-1])
    save_source(output,'commonality',components=parts,coefficients=coef,target=target,predicted=pred)
    gradients=rng.normal(size=(100,6));beta=partial_map_coefficients(target,gradients)
    save_source(output,'gradient_alignment',effect=target,gradients=gradients,beta=beta)
    x=np.linspace(1,3,100);y=.2+1.1*x;estimated=fit_predict(x,y,x,'Affine')
    np.testing.assert_allclose(y,estimated,atol=1e-10)
    save_source(output,'traveling_calibration',observed=y,predicted=estimated)
    cells,counts=grid_overlap(rng.uniform(0,1e6,(100,2)),rng.uniform(0,1e6,(50,2)),100000)
    save_source(output,'geographic_overlap',cells=cells,counts=counts)
    print('ADDITIONAL_METHODS_PASS',output)


if __name__=='__main__':run(Path(__file__).resolve().parents[1]/'outputs/additional_methods')
