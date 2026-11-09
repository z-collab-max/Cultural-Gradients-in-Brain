"""Rebuild study thickness effects and direct province shuffles from brain fields."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from scipy.io import savemat
from matlab_source_reader import read_mat
from china_mixed import CONTROLS, z, read_design, fit_many, province_draws
DIMENSIONS=['Hofstede_indvres_control','Hofstede_powerdis_control','Hofstede_unceravoi_control','z_tightness_control','Hofstede_indvcoll_control','Hofstede_masvfem_control','Hofstede_lonvsho_control']
CONTROL=CONTROLS

def run(brain,culture,out,family,draws,seed,index_file):
 out.mkdir(parents=True,exist_ok=True)
 roi=[n.replace('_area','_thickness') for n in json.loads(Path(__file__).with_name('roi_names_308.json').read_text())['area']]
 cul=read_mat(culture,'Culdata')
 if family=='pc7':fields=[f'pc7_{j}' for j in [1,2,3]];cf=fields
 elif family=='hof6':fields=[f'ValidP_hofPC6_{j}_birth' for j in [1,2,3]];cf=[f'ValidP_hofPC6_{j}' for j in [1,2,3]]
 else:fields=DIMENSIONS;cf=fields
 with h5py.File(brain) as h:
  d=h['brain_data'];cols=lambda keys:np.column_stack([d[k][()].ravel() for k in keys])
  X,groups=read_design(d);Y=cols(roi);B=cols(fields);center=cols(['centerID_rerange']).ravel()
  valid=np.isfinite(np.c_[X,Y,B,center]).all(1)
  X=X[valid];Y=Y[valid];B=B[valid];groups=groups[valid];center=center[valid]
 _,center=np.unique(center,return_inverse=True);sums=np.zeros((center.max()+1,len(roi)));np.add.at(sums,center,Y)
 Y=z(Y-(sums/np.bincount(center)[:,None])[center]);n=len(Y)
 _,first,group=np.unique(groups[:,0],return_index=True,return_inverse=True)
 indices=None
 if index_file:
  fixed=np.load(index_file);indices=fixed['indices_zero_based'][:,:draws].T
  culture_values=np.column_stack([np.asarray(cul[k]).ravel() for k in cf])
  if not np.allclose(culture_values,B[first]):raise ValueError('Culture rows must match sorted province identifiers')
 maps=[];nulls=[]
 for j,field in enumerate(fields):
  predictors,indices=province_draws(B[:,j],groups,draws,seed,indices)
  result=fit_many(X,Y,predictors,groups);beta=result['beta'][1:];tv=result['t'][0]
  data=dict(observed=result['beta'][0],t=tv,p=result['p'][0],df=result['df'][0],null_beta=beta,n=n,roi_names=np.array(roi),predictor=field,shuffle_indices=indices.T)
  np.savez_compressed(out/(field+'.npz'),**data)
  savemat(out/(field+'.mat'),{k:v.astype(object) if isinstance(v,np.ndarray) and v.dtype.kind=='U' else v for k,v in data.items()},do_compression=True)
  maps.append(result['beta'][0]);nulls.append(beta);print(family,field,n,draws,flush=True)
 np.savez_compressed(out/'models.npz',maps=np.array(maps),nulls=np.array(nulls),predictor_names=np.array(fields),roi_names=np.array(roi),controls=np.array(CONTROL))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for key in ['brain','culture','output']:p.add_argument('--'+key,required=True,type=Path)
 p.add_argument('--family',choices=['pc7','hof6','dimensions'],required=True)
 p.add_argument('--indices',type=Path);p.add_argument('--draws',type=int,default=10000);p.add_argument('--seed',type=int,default=42)
 a=p.parse_args();run(a.brain,a.culture,a.output,a.family,a.draws,a.seed,a.indices)
