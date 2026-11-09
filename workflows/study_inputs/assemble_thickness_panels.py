"""Create eigenmode/PC3 panels or Hofstede six-gradient and Webb panels."""
import argparse,sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat,savemat
from matlab_source_reader import read_mat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from panel_data import panel,save
from assemble_alternative_gradients import build,NAMES

def z(x):return (x-x.mean(0))/x.std(0,ddof=1)
def tail(a,b):return (abs(a)>=abs(b)).mean(0)
def run(source,gradients,spins,out,family,webb):
 out.mkdir(parents=True,exist_ok=True)
 with np.load(source/'models.npz') as f:d=dict(f)
 G=np.column_stack([read_mat(gradients,k) for k in ['SAaxis','mode2','reardon_dev','xu_evo','intervar','cogPC1']]);G[:,1]*=-1;G=z(G)
 if family=='pc7':
  eig=z(read_mat(gradients,'DK308_modes_asym')[:,1:20]);b=np.linalg.pinv(eig)@z(d['maps'][0]);nb=z(d['nulls'][0].T).T@np.linalg.pinv(eig).T
  data=dict(modes=np.arange(2,21),reference_maps=eig,effect=z(d['maps'][0]),beta=b,null_beta=nb,p_two_sided=tail(nb,b),p_min_tail=tail(nb,b)/2,display_order_0based=np.argsort(-abs(b)))
  save('ExtendedDataFigure7',out/'eigenmodes',data,[panel('ExtendedDataFigure7','Eigenmodes','vertical_violin',**data)])
  b=np.linalg.pinv(G)@z(d['maps'][2]);nb=z(d['nulls'][2].T).T@np.linalg.pinv(G).T
  with np.load(source/'pc7_3.npz') as f:tv=f['t'];df=f['df']
  data=dict(t=tv,df=df,surface_variance_estimator='REML mixed model',names=np.array(NAMES),partial_beta=b,province_null=nb,p_province=tail(nb,b))
  save('ExtendedDataFigure7',out/'pc3',data,[panel('ExtendedDataFigure7','Culture_PC3','violin',beta=b,null=nb,p=data['p_province'],labels=np.array(NAMES,object))])
  savemat(out/'PC3_surface.mat',dict(roi_names=d['roi_names'].astype(object),values=tv,color_limits=[-4,4]))
 else:
  if not webb:raise ValueError('Hofstede model requires the generated Webb result')
  w=loadmat(webb,simplify_cells=True);idx=read_mat(spins,'perm_id').astype(int)-1
  main=dict(beta=d['maps'][:2],null_beta=d['nulls'][:2].transpose(1,0,2),cluster_influence=w['influence_all'][:2],webb_weights=w['multG_all_cn'],roi_names=d['roi_names'])
  spatial=dict(gradients=G,gradient_labels=NAMES,spin_indices=idx.T,roi_names=d['roi_names'])
  for j,obs in enumerate(w['Obs'][:2]):np.testing.assert_allclose(d['maps'][j]*obs['sidir'],obs['b_aligned'],atol=1e-10,rtol=0)
  build(main,spatial,out,10000,42)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['source','gradients','output']:p.add_argument('--'+k,type=Path,required=True)
 p.add_argument('--spins',type=Path);p.add_argument('--webb',type=Path);p.add_argument('--family',choices=['pc7','hof6'],required=True)
 a=p.parse_args();run(a.source,a.gradients,a.spins,a.output,a.family,a.webb)
