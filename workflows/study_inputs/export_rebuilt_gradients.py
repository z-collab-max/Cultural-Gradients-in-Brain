"""Export ED3 panel data from rebuilt hemisphere-derived gradient inputs."""
from pathlib import Path
import argparse
import numpy as np
from scipy.io import loadmat,savemat
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
a=parser.parse_args();a.output.mkdir(parents=True,exist_ok=True)
d=loadmat(a.inputs,simplify_cells=True)
maps=np.column_stack([d['SAaxis'],-d['DK308_modes_asym'][:,1],d['xu_evo'],d['reardon_dev'],d['intervar'],d['cogPC1'],d['myelin'],d['gene_PC1'],d['fcgradientsM'][:,0]])
names=np.array(['S-A Axis','Eigenmode 2','Evolutionary Exp','Developmental Exp','Intersubject Var','Cognition PC1','Myelin Content','Gene PC1','Functional Grad 1'])
values=dict(maps=maps,names=names,correlations=np.corrcoef(maps.T))
np.savez_compressed(a.output/'EDFig3_gradients.npz',**values)
savemat(a.output/'EDFig3_gradients.mat',{**values,'names':names.astype(object)},do_compression=True)
print('Saved nine ED3 gradient maps and their correlation matrix.')
