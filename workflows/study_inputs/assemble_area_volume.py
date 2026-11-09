"""Assemble 307-ROI area/volume gradient tests with full-atlas rotations."""
import argparse
from pathlib import Path
import sys
import numpy as np
from scipy.io import savemat
from matlab_source_reader import read_mat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from panel_data import save,panel

def z(x):return (x-x.mean(0))/x.std(0,ddof=1)
def tail(null,point):return (abs(null)>=abs(point)).mean(0)
def run(source,gradients,spins,out):
    G=np.column_stack([read_mat(gradients,k) for k in ['SAaxis','mode2','reardon_dev','xu_evo','intervar','cogPC1']]);G[:,1]*=-1
    indices=np.asarray(read_mat(spins,'perm_id'))
    if not np.equal(indices,np.floor(indices)).all() or indices.min()<1 or indices.max()>len(G):raise ValueError('Expected one-based atlas rotations')
    indices=indices.astype(int)-1
    if indices.shape[0]!=len(G):raise ValueError('Rotations must have ROI rows')
    names=np.array(['S-A','Eigenmode2','Developmental','Evolutionary','Intersubject','CognitivePC1'])
    out.mkdir(parents=True,exist_ok=True)
    for modality in ['area','volume']:
        for pc in [1,2]:
            with np.load(source/f'{modality}_PC{pc}.npz') as f:d=dict(f)
            count=len(d['roi_names'])
            if not np.array_equal(d['roi_names'],d['spin_source_roi_names'][:count]):raise ValueError('Selected ROIs must be the ordered atlas prefix')
            sign=1 if pc==1 else -1;full=d['spin_source_effect']*sign;y=full[:count]
            g=z(G[:count]);projector=np.linalg.pinv(g);beta=projector@z(y)
            rotations=z(full[indices[:count]]);sn=(projector@rotations).T
            pn=(projector@z(d['null'].T*sign)).T
            data=dict(roi_names=d['roi_names'],gradients=g,effect=y,t=d['t']*sign,partial_beta=beta,province_null=pn,p_province=tail(pn,beta),spin_null=sn,p_spin=tail(sn,beta),names=names,spin_indices=indices[:count]+1,spin_source_effect=full)
            cols=np.arange(6) if pc==1 else np.array([0,2,3,4,5]);sg=g[:,cols];sb=np.linalg.pinv(sg)@z(y);ss=(np.linalg.pinv(sg)@rotations).T
            data.update(display_spin_beta=sb,display_spin_null=ss,display_spin_names=names[cols],display_spin_p=tail(ss,sb))
            stem=f'{modality}_PC{pc}';fig='ExtendedDataFigure7'
            panels=[panel(fig,stem+'_province','violin',beta=beta,null=pn,p=data['p_province'],labels=names.astype(object)),panel(fig,stem+'_spin','violin',beta=sb,null=ss,p=data['display_spin_p'],labels=names[cols].astype(object))]
            save(fig,out/stem,data,panels)
            savemat(out/(stem+'_surface.mat'),dict(roi_names=d['roi_names'].astype(object),values=data['t'],color_limits=np.array([-1,1])*np.ceil(abs(data['t']).max())))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','gradients','spins','output']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();run(a.source,a.gradients,a.spins,a.output)
