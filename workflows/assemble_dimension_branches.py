"""Assemble dimension-omission and effect-PCA panels."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analyses'))
import numpy as np
from scipy.io import loadmat
from specialized import effect_pca,heat
from map_statistics import spin_test
from panel_data import save

def run(models,spatial,out,components):
 d=loadmat(models,simplify_cells=True)
 with np.load(spatial) as f:s=dict(f)
 if not np.array_equal(np.asarray(d['roi_names']).astype(str),s['roi_names']):raise ValueError('ROI names/order mismatch')
 spin=s['spin_indices']
 if spin.dtype.kind not in 'iu' or spin.shape[1]!=len(s['roi_names']) or spin.min()<0 or spin.max()>=len(s['roi_names']):raise ValueError('Expected zero-based rotation-by-ROI indices')
 labels=np.asarray(d['dimension_labels']).astype(str)
 tests=[[spin_test(d['full_maps'][j],fold[j],spin) for j in range(2)] for fold in d['fold_maps']]
 correlations=np.array([[t[0] for t in row] for row in tests]);pvalues=np.array([[t[2] for t in row] for row in tests])
 source=dict(full_maps=d['full_maps'],fold_maps=d['fold_maps'],correlations=correlations,p_spin=pvalues,dimension_labels=labels)
 panels=[heat('CulturalDimensionOmission','China_dimension_omission',correlations,labels,['PC1','PC2'],pvalues)]
 source.update(matching=d['matching'],orientation=d['orientation'],roi_names=s['roi_names'])
 save('CulturalDimensionOmission',out,source,panels)
 source,panels,_=effect_pca(dict(dimension_effects=d['dimension_effects'],culture_effects=d['full_maps'].T,dimension_labels=labels,spin_indices=spin,reference_component=components,left_hemisphere=np.char.startswith(s['roi_names'],'lh_')))
 source.update(roi_names=s['roi_names'],dimension_effects=d['dimension_effects'],culture_effects=d['full_maps'].T)
 save('ExtendedDataFigure15',out,source,panels)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['models','spatial','output']:p.add_argument('--'+k,type=Path,required=True)
 p.add_argument('--components',type=int,nargs=2,required=True)
 a=p.parse_args();run(a.models,a.spatial,a.output,a.components)



