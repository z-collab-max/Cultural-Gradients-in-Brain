"""Assemble China dimension correspondence and standardized six-map PCA."""
import argparse,sys
from pathlib import Path
import numpy as np
from scipy.io import savemat
from matlab_source_reader import read_mat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analyses'))
from panel_data import save,panel
from specialized import heat
from cultural_pca import cultural_pca

def unit(x):
 x=x-x.mean(0);return x/np.linalg.norm(x,axis=0)
def z(x):return (x-x.mean(0))/x.std(0,ddof=1)
def run(dimensions,hofstede,gradients,spins,out):
 with np.load(dimensions/'models.npz') as f:d=dict(f)
 with np.load(hofstede/'models.npz') as f:h=dict(f)
 if not np.array_equal(d['roi_names'],h['roi_names']):raise ValueError('ROI names/order mismatch')
 G=np.column_stack([read_mat(gradients,k) for k in ['mode2','SAaxis','reardon_dev','xu_evo','intervar','cogPC1']]);G[:,0]*=-1
 maps=np.c_[d['maps'].T,G];labels=np.array(['IVR','PDI','UAI','TL','IDV','MAS','LTO','eigenmode2','SAaxis','devexp','evoexp','intersubjvar','cogpc1'])
 spin=np.asarray(read_mat(spins,'perm_id'),int)-1;Z=unit(maps);R=Z.T@Z;ps=np.ones((13,13));pp=np.full((13,13),np.nan)
 spin_null=np.empty((spin.shape[1],13,13));permutation_null=np.full((d['nulls'].shape[1],13,13),np.nan)
 for j in range(13):
  null=unit(maps[:,j][spin]).T@Z;spin_null[:,j,:]=null;ps[j]=(abs(null)>=abs(R[j])).mean(0)
 for j in range(7):
  null=unit(d['nulls'][j].T).T@Z;permutation_null[:,:,j]=null;pp[:,j]=(abs(null)>=abs(R[:,j])).mean(0)
 for i in range(7):
  for j in range(7,13):pp[i,j]=pp[j,i];permutation_null[:,i,j]=permutation_null[:,j,i]
 pf=np.where(np.isnan(pp),ps,np.maximum(pp,ps));np.fill_diagonal(pf,np.nan)
 recipe=heat('ExtendedDataFigure13','China_dimension_gradients',R,labels,p=pf)
 recipe.update(spin_null=spin_null,permutation_null=permutation_null,symmetric_spin=False,symmetric_permutation=False,spin_only_cells=True)
 save('ExtendedDataFigure13',out,dict(spin_null=spin_null,permutation_null=permutation_null,R=R,p_spin=ps,p_province=pp,p_combined=pf,labels=labels,maps=maps,roi_names=d['roi_names']),[recipe])
 chosen=[0,1,2,4,5,6];raw=d['maps'][chosen].T;fit=cultural_pca(raw)
 for j in range(2):
  if np.corrcoef(fit['scores'][:,j],raw[:,0])[0,1]<0:fit['scores'][:,j]*=-1;fit['loadings'][:,j]*=-1
 culture=h['maps'][:2].T;signs=np.array([np.sign(np.corrcoef(fit['scores'][:,j],culture[:,j])[0,1]) for j in range(2)]);display=culture*signs
 data=dict(scores=fit['scores'],loadings=fit['loadings'],explained=fit['explained']*100,labels=labels[chosen],culture_maps=culture,culture_display=display,display_sign=signs,standardized_maps=z(raw),raw_maps=raw,roi=d['roi_names'])
 panels=[panel('ExtendedDataFigure15',f'China_PC{j+1}_effect','identity_scatter',x=z(display[:,j]),y=z(fit['scores'][:,j]),xlabel=f'Culture PC{j+1} effect',ylabel=f'Effect PC{j+1}',stats='r = %.3f'%np.corrcoef(display[:,j],fit['scores'][:,j])[0,1]) for j in range(2)]
 panels.append(panel('ExtendedDataFigure15','China_effect_biplot','pca_biplot',scores=fit['scores'][:,:2],loadings=fit['loadings'][:,:2],labels=labels[chosen].astype(object),explained=data['explained'][:2]))
 panels.append(panel('ExtendedDataFigure15','China_effect_scree','bars',y=data['explained'],labels=np.array([f'PC{j+1}' for j in range(6)],object),ylabel='Variance (%)'))
 save('ExtendedDataFigure15',out,data,panels)
 for j in range(2):savemat(out/f'Effect_PC{j+1}_surface.mat',dict(roi_names=d['roi_names'].astype(object),values=z(fit['scores'][:,j]),color_limits=[-4,4]))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['dimensions','hofstede','gradients','spins','output']:p.add_argument('--'+k,type=Path,required=True)
 a=p.parse_args();run(a.dimensions,a.hofstede,a.gradients,a.spins,a.output)
