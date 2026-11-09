"""Aggregate participant origins and sites for the three-cohort geography panels."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from geography_inputs import load_china_data,load_usa
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'figures'))
from render_maps import render


def run(a):
    a.output.mkdir(parents=True,exist_ok=True)
    china=load_china_data(a.brain,ROOT/'assets/geography/china_shapes/bou2_4p_utf8.shp');china.pop('province_shapes')
    participants=pd.read_csv(a.nations_participants,dtype={'row_key':str,'unit':str,'batch':str})
    if participants.row_key.duplicated().any():raise ValueError('Duplicated participant row key')
    geo=pd.read_csv(a.nations_coordinates).rename(columns={'nation':'unit','cap_lat':'latitude','cap_lon':'longitude'})
    geo['unit']=geo.unit.replace({'South Korea':'SouthKorea','Czech Republic':'Czechia'})
    if geo.unit.duplicated().any():raise ValueError('Duplicated unit coordinates')
    geo=geo.set_index('unit');nations=[]
    for unit,rows in participants.groupby('unit'):
        loc=geo.loc[unit];nations.append(dict(name=unit,n=len(rows),sites=rows.batch.nunique(),lat=float(loc.latitude),lon=float(loc.longitude)))
    usa=load_usa(a.usa_participants)
    data=dict(china={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in china.items()},nations=nations,usa=usa)
    (a.output/'source_data.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    centers=[dict(longitude=float(lon),latitude=float(lat),count=1,marker='*') for lon,lat in zip(china['center_lon'],china['center_lat'])]
    origins=[dict(longitude=float(lon),latitude=float(lat),count=int(count),composition=row.tolist()) for lon,lat,count,row in zip(china['prov_lon'],china['prov_lat'],china['totals'],china['counts'])]
    panels=[dict(region='china',title='Scanning centers',extent=[73,135,15,54],points=centers)]
    for region,rows in [('world',nations),('usa',usa)]:panels.append(dict(region=region,title='Scanning centers',points=[dict(longitude=r['lon'],latitude=r['lat'],count=r['sites'],marker='*') for r in rows]))
    panels.append(dict(region='china',title='Participant origins',extent=[73,135,15,54],points=origins))
    for region,rows in [('world',nations),('usa',usa)]:panels.append(dict(region=region,title='Participant origins',points=[dict(longitude=r['lon'],latitude=r['lat'],count=r['n']) for r in rows]))
    config=dict(columns=3,panels=panels);(a.output/'maps.json').write_text(json.dumps(config,indent=2));render(config,a.output/'center_geography')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['brain','nations-participants','nations-coordinates','usa-participants','output']:p.add_argument('--'+key,type=Path,required=True)
    run(p.parse_args())
