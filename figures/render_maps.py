"""Render aggregated geographic panels using bundled boundary geometry."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from map_geometry import albers,robinson
ASSETS=Path(__file__).resolve().parents[1]/'assets/geography'

def draw(ax,panel):
    region=panel['region']
    if region=='china':files=['china_nat_line.json','china_prov_line.json'];project=lambda x,y:albers(x,y,105,0,25,47)
    elif region=='usa':files=['usa.json'];project=lambda x,y:albers(x,y,-96,23,29.5,45.5)
    elif region=='world':files=['world.json'];project=robinson
    else:raise ValueError('region must be china, usa or world')
    values=panel.get('region_values',{});cmap=plt.get_cmap('RdBu_r');limits=panel.get('color_limits',[-1,1]);norm=plt.Normalize(*limits)
    if 'colors' in panel:
        from matplotlib.colors import LinearSegmentedColormap
        cmap=LinearSegmentedColormap.from_list('region_scores',panel['colors'])
    if region=='china' and values:raise ValueError('China boundary-line panels use point values; polygon fill requires a named polygon atlas')
    known=set()
    for filename in files:
        geo=json.loads((ASSETS/filename).read_text(encoding='utf-8'))
        for feature in geo['features']:
            g=feature['geometry'];name=feature['properties'].get('name','');color=cmap(norm(values[name])) if name in values else '#eeeeee'
            known.add(name)
            if name in panel.get('exclude_regions',[]) or ('only_regions' in panel and name not in panel['only_regions']):continue
            if g['type']=='MultiPolygon':polys=g['coordinates']
            elif g['type']=='Polygon':polys=[g['coordinates']]
            elif g['type']=='MultiLineString':polys=[[v] for v in g['coordinates']]
            else:polys=[[g['coordinates']]]
            for poly in polys:
                for k,ring in enumerate(poly):
                    xy=np.asarray(ring,float)
                    if len(xy)<2:continue
                    x,y=project(xy[:,0],xy[:,1])
                    if g['type'].endswith('Polygon'):ax.add_patch(Polygon(np.c_[x,y],facecolor=color if k==0 else 'white',edgecolor='#aaaaaa',linewidth=.35))
                    else:ax.plot(x,y,color='#aaaaaa',linewidth=.35)
            if name in panel.get('abbreviations',{}) and g['type'].endswith('Polygon'):
                from shapely.geometry import shape
                geometry=shape(g);parts=list(geometry.geoms) if geometry.geom_type=='MultiPolygon' else [geometry]
                label=max(parts,key=lambda part:part.area).representative_point();px,py=project(label.x,label.y)
                ax.text(px,py,panel['abbreviations'][name],fontsize=5,ha='center',va='center')
    if set(values)-known:raise ValueError('Region names not found: '+str(sorted(set(values)-known)))
    if region in ['china','world']:
        for line in json.loads((ASSETS/'nine_dash.json').read_text(encoding='utf-8')):
            xy=np.asarray(line);x,y=project(xy[:,0],xy[:,1]);ax.plot(x,y,color='#aaaaaa',linewidth=.35)
    for row in panel.get('points',[]):
        x,y=project(row['longitude'],row['latitude']);size=10+5*np.sqrt(row.get('count',1))
        if 'composition' in row:
            a=np.asarray(row['composition'],float)
            if (a<0).any() or a.sum()<=0:raise ValueError('Invalid composition')
            start=0
            for j,share in enumerate(a/a.sum()):
                t=np.linspace(start,start+2*np.pi*share,25);marker=np.vstack([[0,0],np.c_[np.cos(t),np.sin(t)],[0,0]])
                ax.scatter([x],[y],s=size,marker=marker,color=plt.get_cmap('tab20')(j%20),linewidths=.2,edgecolors='white');start+=2*np.pi*share
        else:ax.scatter([x],[y],s=size,marker=row.get('marker','o'),color=row.get('color','#71819e'),edgecolors='white',linewidths=.5)
    if 'extent' in panel:
        lo,hi,bottom,top=panel['extent'];xx,yy=project(np.array([lo,lo,hi,hi]),np.array([bottom,top,bottom,top]));ax.set_xlim(min(xx),max(xx));ax.set_ylim(min(yy),max(yy))
    else:ax.autoscale_view()
    ax.set_aspect('equal');ax.axis('off');ax.set_title(panel.get('title',''),fontsize=8)
    if values and panel.get('colorbar',True):plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,shrink=.55)
    if 'legend' in panel:
        handles=[ax.scatter([],[],s=v.get('size',25),color=v['color'],label=v['label']) for v in panel['legend']]
        ax.legend(handles=handles,loc=panel.get('legend_location','lower left'),frameon=False)
    for child in panel.get('insets',[]):
        inset=ax.inset_axes(child['position']);draw(inset,child['panel'])

def render(config,output):
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':.85})
    panels=config['panels'];cols=int(config.get('columns',min(3,len(panels))));rows=(len(panels)+cols-1)//cols
    fig,axes=plt.subplots(rows,cols,figsize=(5*cols,3.5*rows),squeeze=False)
    for ax,p in zip(axes.flat,panels):draw(ax,p)
    for ax in list(axes.flat)[len(panels):]:ax.axis('off')
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True);fig.tight_layout()
    fig.savefig(output.with_suffix('.svg'));fig.savefig(output.with_suffix('.png'),dpi=180);plt.close(fig)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(json.loads(a.input.read_text(encoding='utf-8')),a.output)
