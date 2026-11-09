"""Miller-grid overlap and clipped geographic source data."""
from pathlib import Path
import argparse,json,sys
import numpy as np
from scipy.io import savemat
from shapely.geometry import Polygon,Point,box,mapping
from shapely.ops import unary_union
from shapely import make_valid
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'figures'))

def miller(lonlat):
    a=np.asarray(lonlat,float);W=6381372*np.pi*2/1000;H=W/2
    y1=1.25*np.log(np.tan(np.pi/4+.4*np.deg2rad(a[:,1])))
    return np.c_[a[:,0]*W/360+W/2,H/2-H*y1/(2*2.3)]

def inverse(xy):
    a=np.asarray(xy,float);W=6381372*np.pi*2/1000;H=W/2;y1=(H/2-a[:,1])*2*2.3/H
    return np.c_[(a[:,0]-W/2)*360/W,np.rad2deg((np.arctan(np.exp(y1/1.25))-np.pi/4)/.4)]

def generate(imaging,culture,size=100):
    if not np.isfinite(size) or size<=0:raise ValueError('Positive grid size required')
    geo=json.loads((ROOT/'assets/geography/china_prov_poly.json').read_text(encoding='utf-8'))
    region=unary_union([make_valid(Polygon(r)) for f in geo['features'] for r in f['geometry']['coordinates'] if len(r)>=4])
    origin=miller([[80,60]])[0];upper=miller([[140,5]])[0];end=origin+np.ceil((upper-origin)/size)*size
    inputs=[np.asarray(imaging,float),np.asarray(culture,float)];counts=[];removed=[]
    for points in inputs:
        if points.ndim!=2 or points.shape[1]!=2 or not np.isfinite(points).all():raise ValueError('Finite longitude/latitude pairs required')
        good=np.array([region.covers(Point(p)) for p in points]);xy=miller(points);good&=((xy>=origin)&(xy<=end)).all(1)
        cells=np.minimum(np.floor((xy[good]-origin)/size).astype(int),np.ceil((upper-origin)/size).astype(int)-1)
        keys,n=np.unique(cells,axis=0,return_counts=True);counts.append({tuple(k):int(v) for k,v in zip(keys,n)});removed.append(int((~good).sum()))
    keys=sorted(set(counts[0])|set(counts[1]));features=[];values=[];category_counts=np.zeros(3,int)
    for key in keys:
        a,b=[c.get(key,0) for c in counts];category=0 if a and b else 1 if a else 2;category_counts[category]+=1
        lo=origin+np.asarray(key)*size;hi=lo+size;corners=inverse(np.array([lo,[hi[0],lo[1]],hi,[lo[0],hi[1]]]))
        clipped=make_valid(Polygon(corners)).intersection(region)
        features.append(dict(type='Feature',geometry=mapping(clipped),properties=dict(category=int(category),imaging=a,culture=b)));values.append([a,b])
    source=dict(grid_cells=np.array(keys,dtype=int).reshape(-1,2),counts=np.array(values,dtype=int).reshape(-1,2),category_counts=category_counts,excluded_counts=removed,grid_size_km=size,miller_origin=origin)
    return source,dict(type='FeatureCollection',features=features)

def render(geo,source,stem):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as Patch
    from render_maps import draw
    from map_geometry import albers
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':.85});fig,ax=plt.subplots(figsize=(8,7));draw(ax,dict(region='china'))
    colors=['#95d1a9','#e5d9b3','#e2bdb7']
    for feature in geo['features']:
        g=feature['geometry'];polys=g['coordinates'] if g['type']=='MultiPolygon' else [g['coordinates']] if g['type']=='Polygon' else []
        for p in polys:
            xy=np.array(p[0]);x,y=albers(xy[:,0],xy[:,1],105,0,25,47);ax.add_patch(Patch(np.c_[x,y],facecolor=colors[feature['properties']['category']],edgecolor='#888888',linewidth=.4,alpha=.92))
    handles=[Patch([[0,0]],facecolor=c,label=f'{label} (N={n})') for c,label,n in zip(colors,['Both','Imaging only','Culture only'],source['category_counts'])]
    ax.legend(handles=handles,loc='lower left',frameon=False);ax.set_title(f"Spatial overlap ({source['grid_size_km']:g} km grid)")
    fig.savefig(str(stem)+'.svg',bbox_inches='tight');fig.savefig(str(stem)+'.png',dpi=180,bbox_inches='tight');plt.close(fig)

def run(imaging,culture,output,size=100):
    output=Path(output);output.mkdir(parents=True,exist_ok=True);s,g=generate(imaging,culture,size)
    np.savez_compressed(output/'source_data.npz',**s);savemat(output/'source_data.mat',s);(output/'grid_polygons.json').write_text(json.dumps(g),encoding='utf-8');render(g,s,output/'grid_overlap');return s
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--size',type=float,default=100);a=p.parse_args()
    with np.load(a.input,allow_pickle=False) as d:run(d['imaging_lonlat'],d['culture_lonlat'],a.output,a.size)
