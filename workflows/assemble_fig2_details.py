"""Assemble Fig. 2 target-ROI scatter, province SNR and replacement panels."""
import argparse
import csv
import json
import os
from pathlib import Path

if __name__=='__main__':
    for setting in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:
        os.environ[setting]='2'
import numpy as np
from scipy.io import loadmat
from fig2_province_model import recovery_draws
from fig2_replacement import replacement
from panel_data import panel, save
from fig2_mixed_adjustment import adjusted_values

TARGETS = ['lh_superiorfrontal_part8_thickness','lh_lateraloccipital_part7_thickness',
           'lh_superiorfrontal_part9_thickness','rh_inferiorparietal_part4_thickness']


def z(x, ddof=0):
    sd=x.std(0,ddof=ddof)
    if np.any(sd<=0):raise ValueError('Constant input to standardization')
    return (x-x.mean(0))/sd


def center_adjust(raw,scanner):
    adjusted=raw.copy()
    for center in np.unique(scanner):
        rows=scanner==center
        adjusted[rows]-=adjusted[rows].mean(0)
    return adjusted


def scatter(raw,C,pcs,province,city,scanner,names,targets,annotations):
    selected=[names.index(name) for name in targets]
    adjusted=adjusted_values(raw[:,selected],C,pcs,province,city,scanner,[1,1,2,2]);points=[];panels=[]
    for j,(key,name) in enumerate(zip(['pc1_pos','pc1_neg','pc2_pos','pc2_neg'],targets)):
        ix=names.index(name);pc=j//2;x=z(pcs[:,pc])*(1 if pc==0 else -1)
        v=adjusted[:,j]
        rows=[]
        for group in np.unique(province):
            hit=province==group
            if np.ptp(x[hit])>1e-10:raise ValueError('Culture must be constant within province')
            rows.append((x[hit][0],v[hit].mean(),v[hit].std()/np.sqrt(hit.sum()),str(group),int(hit.sum())))
        rows.sort();xs,means,se,groups,ns=zip(*rows)
        points.extend((key,name,i+1,*row[:3],len(x),row[3],row[4]) for i,row in enumerate(rows))
        beta,p=annotations[j]
        panels.append(panel('Figure2','Fig2b_'+key,'error_scatter',x=np.array(xs),y=np.array(means),se=np.array(se),
            color=[.72,.50,.46] if j%2==0 else [.44,.51,.62],stats=fr'\beta = {beta:.3g}, p = {p:.3g}',
            xlabel=f'Culture PC{pc+1} (z)',ylabel='Adjusted thickness (mm)'))
    return points,panels


def calibrate(cal,val,sd,beta,se,sigma,tau,x):
    scale=se*sd/cal.std(0,ddof=1)
    if not np.isfinite(scale).all() or np.any(scale<=0):raise ValueError('Invalid SE calibration')
    signal=x.var()*np.mean(beta**2);noise=np.mean(scale**2*(sigma+tau)/sd**2)
    return dict(truth=beta*sd,source_se=se*sd,sd_um=sd,beta_null=val*scale,scale=scale,
                baseline=signal/noise,signal_variance=signal,noise_variance=noise)


def snr(samples):
    x=np.unique(np.r_[np.geomspace(1e-5,1e-1,201),[d['baseline'] for d in samples]])
    out={}
    for k,d in enumerate(samples,1):
        t=d['truth']/d['sd_um'];e=d['beta_null']/d['sd_um'];t=t-t.mean();e=e-e.mean(1,keepdims=True)
        c=np.sqrt(d['baseline']/x);cross=e@t;s=t@t;en=np.sum(e*e,axis=1)
        r=(s+c[:,None]*cross)/np.sqrt(s*(s+2*c[:,None]*cross+c[:,None]**2*en))
        q=np.quantile(r*r,[.025,.5,.975],axis=1);K=np.mean(t*t)/np.mean(np.var(e,axis=0,ddof=1))/d['baseline'];tx=np.geomspace(x[0],x[-1],500)
        vals=dict(x=x,lower=q[0],median=q[1],upper=q[2],theory_x=tx,theory_y=K*tx/(1+K*tx),actual_x=d['baseline'],actual_q=q[:,np.flatnonzero(x==d['baseline'])[0]],noise_scales=c)
        out.update({f'pc{k}_{name}':v for name,v in vals.items()})
    return out


def write_sources(output,points,panels,samples,test_orders):
    dest=save('Figure2',output/'roi',dict(annotations=np.array([p['stats'] for p in panels]),x=[v[3] for v in points],mean_mm=[v[4] for v in points],sem_mm=[v[5] for v in points],
        panel=np.array([v[0] for v in points]),roi=np.array([v[1] for v in points]),province=np.array([v[7] for v in points]),n=[v[8] for v in points]),panels)
    with (dest/'province_points.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['panel','ROI','order','culture_PC_z','mean_mm','SEM_mm','total_N','province','province_n']);writer.writerows(points)
    (dest/'render_manifest.json').write_text(json.dumps([dict(kind='native',input='panels.mat')]))
    ablation,tests=replacement(samples,test_orders)
    for key,values,layout in [('snr',snr(samples),'plot_fig2e_snr'),('replacement',ablation,'plot_fig2e_replacement')]:
        dest=save('Figure2',output/key,values,[])
        (dest/'render_manifest.json').write_text(json.dumps([dict(kind='layout',input='source_data.mat',function=layout)]))
    (dest/'statistics.json').write_text(json.dumps(tests,indent=2,allow_nan=False))


def build(results,output,targets,draws,test_orders,scatter_results=None,annotation_file=None):
    if draws<3 or test_orders<2:raise ValueError('Use at least 3 draws per calibration batch and 2 reference orders')
    if len(targets)!=4:raise ValueError('Expected two ROI names for each cultural PC')
    d=loadmat(results/'fig2_individual_inputs.mat',simplify_cells=True);m=loadmat(results/'main_maps.mat',simplify_cells=True)
    if 'city' not in d:raise ValueError('Recovery inputs require birthplace city IDs; regenerate the shared China chain outputs')
    names=np.asarray(d['roi_names']).ravel().tolist()
    if not np.array_equal(names,np.asarray(m['roi_names']).ravel()):raise ValueError('Mismatch in ROI order')
    raw=d['raw_mm'];C=d['covariates'];pcs=d['birth_scores'];province=d['province'];scanner=d['scanner']
    scatter_data,scatter_model=d,m
    if scatter_results is not None:
        scatter_data=loadmat(scatter_results/'fig2_individual_inputs.mat',simplify_cells=True)
        scatter_model=loadmat(scatter_results/'main_maps.mat',simplify_cells=True)
    snames=np.asarray(scatter_data['roi_names']).ravel().tolist()
    if not np.array_equal(snames,np.asarray(scatter_model['roi_names']).ravel()):raise ValueError('Scatter ROI order mismatch')
    annotations=[(scatter_model['beta'][j//2,snames.index(v)],scatter_model['p_model'][j//2,snames.index(v)]) for j,v in enumerate(targets)]
    if annotation_file:
        rows=json.loads(annotation_file.read_text(encoding='utf-8-sig'))
        if [r['panel'] for r in rows]!=['pc1_pos','pc1_neg','pc2_pos','pc2_neg']:raise ValueError('Unexpected annotation order')
        annotations=[(r['beta'],r['p']) for r in rows]
    if 'city' not in scatter_data:raise ValueError('Scatter inputs require birthplace city IDs; regenerate the shared China chain outputs')
    points,panels=scatter(scatter_data['raw_mm'],scatter_data['covariates'],scatter_data['birth_scores'],scatter_data['province'],scatter_data['city'],scatter_data['scanner'],snames,targets,annotations)
    Y=center_adjust(raw*1000,scanner);sd=Y.std(0,ddof=1);C=z(C,1);pcs=z(pcs,1)
    recovery=recovery_draws(Y,C,pcs,province,d['city'],draws)
    samples=[];folder=output/'draws';folder.mkdir(parents=True,exist_ok=True)
    for k in range(2):
        x=pcs[:,k];cal=recovery['calibration_errors'][:,:,k];val=recovery['validation_errors'][:,:,k]
        sigma=recovery['residual_variance'][k];vp=recovery['province_variance'][k];vc=recovery['city_variance'][k]
        sample=calibrate(cal,val,sd,m['beta'][k],m['se'][k],sigma,vp+vc,x);samples.append(sample)
        np.savez_compressed(folder/f'pc{k+1}.npz',**sample,calibration_errors=cal,validation_errors=val,
            residual_variance=sigma,province_variance=vp,city_variance=vc,
            recovery_beta=recovery['observed_beta'][k],recovery_se=recovery['observed_se'][k],
            calibration_seed=recovery['calibration_seed'].ravel()[k],validation_seed=recovery['validation_seed'].ravel()[k])
    write_sources(output,points,panels,samples,test_orders)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--targets',nargs=4,default=TARGETS);p.add_argument('--draws',type=int,default=1000);p.add_argument('--test-orders',type=int,default=2000)
    p.add_argument('--scatter-results',type=Path);p.add_argument('--annotations',type=Path)
    a=p.parse_args();build(a.results,a.output,a.targets,a.draws,a.test_orders,a.scatter_results,a.annotations)
