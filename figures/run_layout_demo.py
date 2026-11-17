"""Generate synthetic inputs for the dedicated layouts, surfaces and maps."""
from pathlib import Path
import argparse,sys,json
import numpy as np
from scipy.io import loadmat,savemat
root=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(root/'workflows'),str(root/'analyses')]
from assemble import assemble
from render_maps import render
from cultural_pca import joint_cultural_pca

def main(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(31)
    def save(name,d):savemat(output/(name+'.mat'),d,do_compression=True)
    names=[]
    for h,prefix in [('L','lh'),('R','rh')]:
        a=loadmat(root/f'assets/surfaces/atlas_{h}.mat',simplify_cells=True)
        names.extend(prefix+'_'+str(n)+'_thickness' for i,n in zip(a['label_ids'],a['label_names']) if i not in [0,8])
    save('brain',dict(roi_names=np.array(names,object),values=np.sin(np.linspace(-3,3,len(names))),color_limits=[-1,1],significant=np.arange(len(names))%50==0))
    save('vertex',dict(values_lh=np.sin(np.linspace(0,6,10242)),values_rh=np.cos(np.linspace(0,6,10242)),color_limits=[-1.,1.]))
    x=rng.normal(size=(30,3));y=rng.normal(size=30)*.1;fit=np.empty((3,100,4))
    for j in range(3):
        xx=np.linspace(-2,2,100);b=np.polyfit(x[:,j],y,1);yy=np.polyval(b,xx);fit[j]=np.c_[xx,yy,yy-.05,yy+.05]
    save('si1',dict(dictionary=np.array([[f'Word{i}_{j}' for j in range(4)] for i in range(10)],object),word_labels=np.array(['W'+str(i) for i in range(20)],object),explained=np.tile(np.arange(10,0,-1)/55*100,(2,1)),fit=fit,x=x,y=y,province_names=np.array(['P'+str(i) for i in range(30)],object),r=np.array([np.corrcoef(v,y)[0,1] for v in x.T]),p=np.ones(3)))
    d=dict(participant1=2.5+rng.normal(size=(10,32))*.1,participant2=2.5+rng.normal(size=(10,32))*.1,center_labels=np.array(['Site'+str(i) for i in range(10)],object))
    s,_=assemble('SupplementaryFigure2',d);save('si2',s)
    d=dict(duration=rng.integers(0,12,100),birth_province=np.array(['P'+str(i%6) for i in range(100)]),destination_province=np.array(['P'+str(i%5) for i in range(100)]))
    s,_=assemble('SupplementaryFigure4',d);s['route_labels']=s['route_labels'].astype(object);save('si4',s)
    snr={}
    for j in [1,2]:
        x=np.logspace(-2,2,20);y=x/(x+1)
        for k,v in dict(x=x,lower=np.maximum(0,y-.05),upper=np.minimum(1,y+.05),median=y,theory_x=x,theory_y=y,actual_q=[.4,.5,.6],actual_x=1).items():snr[f'pc{j}_{k}']=v
    save('snr',snr)
    joint=joint_cultural_pca(rng.normal(size=(30,6)),rng.normal(size=(26,6)),np.zeros(6),np.ones(6)*.2)
    a=joint['score_joint_provinces'];b=joint['score_joint_nations'];coeff=joint['coeff_joint'];explained=joint['explained_joint']
    save('ed15a',dict(scores_china=a,scores_world=b,coeff_joint=coeff,explained_joint=explained,world_nations_clean=np.array(['Nation'+str(i) for i in range(26)],object),prov_display=np.array(['P'+str(i) for i in range(30)],object),continent_grp=np.arange(26)%4+1))
    save('ed15b',joint)
    panels=[]
    for kind in ['Centers','Origins']:
        for region,lon,lat in [('china',110,35),('world',5,45),('usa',-95,38)]:
            points=[dict(longitude=lon+i,latitude=lat+i/2,count=10+i*5) for i in range(3)]
            if kind=='Origins' and region=='china':
                for p in points:p['composition']=[2,3,5]
            panel=dict(region=region,title=kind,points=points)
            if region=='usa':panel['extent']=[-125,-66,24,50]
            panels.append(panel)
    config=dict(columns=3,panels=panels);(output/'maps.json').write_text(json.dumps(config),encoding='utf-8');render(config,output/'maps')
    print(output)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=root/'outputs/layout_demo');a=p.parse_args();main(a.output)
