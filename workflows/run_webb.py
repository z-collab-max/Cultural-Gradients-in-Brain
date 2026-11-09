from pathlib import Path
import argparse,sys
import numpy as np
from scipy.io import savemat
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analyses'));sys.path.insert(0,str(ROOT/'international'))
from webb_statistics import cluster_score_intervals

def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['cluster','international'],required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--repetitions',type=int,default=10000);p.add_argument('--seed',type=int,default=42);a=p.parse_args()
    with np.load(a.input,allow_pickle=False) as f:d=dict(f)
    if a.mode=='cluster':out=cluster_score_intervals(d['point'],d['influence'],d['gradients'],a.repetitions,a.seed,weights=d.get('weights'))
    else:
        from webb_intervals import webb_intervals
        out=webb_intervals(d,a.repetitions,a.seed)
    a.output.mkdir(parents=True,exist_ok=True);np.savez_compressed(a.output/'webb_statistics.npz',**out);savemat(a.output/'webb_statistics.mat',out)
    print(a.output)
if __name__=='__main__':main()
