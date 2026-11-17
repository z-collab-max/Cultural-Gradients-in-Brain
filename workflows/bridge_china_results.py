"""Connect the China model result file to best-order hierarchy source data."""
from pathlib import Path
import argparse
import numpy as np
from scipy.io import loadmat
from run_specialized import run,read_input

def bridge(results,spatial,output,repetitions=10000):
    s=loadmat(results,simplify_cells=True);m=s['mediation'];order=s['behavior_order']
    names=np.empty(len(order['best_order']),object);names[np.asarray(order['best_order'],int)-1]=np.asarray(order['task_names']).ravel()
    effect=np.asarray(m['task_indirect_z'])
    if not np.isfinite(effect).all():effect=np.asarray(m['task_indirect'])
    d=dict(spatial,task_effects=np.stack([m['behavior_effect'],effect[:,:,0],effect[:,:,1]]),task_labels=names,branch_labels=np.array(['Behavior','CulturePC1','CulturePC2']))
    return run('hierarchy',d,output,repetitions)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--spatial',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();bridge(a.results,read_input(a.spatial),a.output)
