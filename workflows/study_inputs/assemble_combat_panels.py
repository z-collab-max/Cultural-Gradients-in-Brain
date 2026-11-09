"""Create separate ComBat source tables, spatial tests and plotting inputs."""
import argparse,sys
from pathlib import Path
import h5py
import numpy as np
from scipy.io import loadmat,savemat
from matlab_source_reader import read_mat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from panel_data import panel,save


def z(x):return (x-x.mean(0))/x.std(0,ddof=1)
def tail(null,b):return (abs(null)>=abs(b)).mean(0)


def run(source,gradients,spins,out,protected):
    m=loadmat(source/f'EDFig5_ComBat_PC{12 if protected==2 else 123}_maps.mat',simplify_cells=True)
    if int(m['protectedPCs'])!=protected:raise ValueError('Observed ComBat protection differs from requested version')
    roi=np.asarray(m['brain_var'],str);names=np.array(['SAaxis','Eigenmode2','Developmental','Evolutionary','Intersubject','CognitivePC1'])
    G=np.column_stack([read_mat(gradients,k) for k in ['SAaxis','mode2','reardon_dev','xu_evo','intervar','cogPC1']]);G[:,1]*=-1;G=z(G)
    spin=read_mat(spins,'perm_id').astype(int)-1;projector=np.linalg.pinv(G)
    with h5py.File(source/'Permutation_Results_FullStats_CulturePC.h5') as h:
        if int(h.attrs['protected_birth_pc_count'])!=protected:raise ValueError('Observed and null ComBat protection must match')
        for pc,sign,primary in [(0,1,1),(1,-1,0)]:
            y=sign*m['coeff'][:,pc];ex=sign*m['expcs'][:,pc];b=projector@z(y)
            nb=z(sign*h['AllPermRes/CombatMod1/Beta'][pc].T).T@projector.T
            keep=np.arange(6) if pc==0 else np.array([0,2,3,4,5]);ps=np.linalg.pinv(G[:,keep]);sb=ps@z(y);sn=z(y[spin]).T@ps.T
            pp=tail(nb,b);sp=tail(sn,sb);pidx=int(np.flatnonzero(keep==primary)[0])
            data=dict(combat=y,expcs=ex,map_r=np.corrcoef(y,ex)[0,1],t=sign*m['tvalue'][:,pc],roi_p=m['pvalue'][:,pc],roi_names=roi,
                gradients=G,names=names,partial_beta=b,province_null=nb,p_province=pp,spin_beta=sb,spin_null=sn,p_spin=sp,spin_names=names[keep],
                protected_birth_pc_count=protected,null_protected_birth_pc_count=int(h.attrs['protected_birth_pc_count']),sample_n=m['sample_n'])
            prefix=f'PC{pc+1}'
            panels=[panel('ExtendedDataFigure5',prefix+'_map_comparison','combat_comparison',x=z(ex),y=z(y),xlabel='Center-adjusted effect',ylabel='ComBat effect'),
                    panel('ExtendedDataFigure5',prefix+'_gradient','scatter',x=G[:,primary],y=z(y),lh=np.char.startswith(roi,'lh_'),xlabel=names[primary],ylabel='Culture effect',stat_null=sn[:,pidx],stat_observed=sb[pidx],stats=f'β = {b[primary]:.3f}; p = {sp[pidx]:.3f}'),
                    panel('ExtendedDataFigure5',prefix+'_province','violin',beta=b,null=nb,p=pp,labels=names.astype(object)),
                    panel('ExtendedDataFigure5',prefix+'_spin','violin',beta=sb,null=sn,p=sp,labels=names[keep].astype(object))]
            folder=save('ExtendedDataFigure5',out/prefix,data,panels)
            savemat(folder/'surface.mat',dict(roi_names=roi.astype(object),values=data['t'],color_limits=[-4,4]))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','gradients','spins','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--protected-pcs',type=int,choices=[2,3],default=2)
    a=p.parse_args();run(a.source,a.gradients,a.spins,a.output,a.protected_pcs)
