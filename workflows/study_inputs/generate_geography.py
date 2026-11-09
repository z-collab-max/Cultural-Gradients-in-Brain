"""Aggregate province/city counts and reconstruct occupied spatial grids."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from primary_residence import extract_primary_residence
from miller_grid import lon_lat_2dikaer,xy_to_coor,get_matlab_grid_idx

def city_counts(culture):
    city=culture.city.astype(str).copy()
    municipal=culture.province.isin(['北京市','天津市','上海市','重庆市'])
    city.loc[municipal]=culture.loc[municipal,'province']
    city=city.replace({'盘州市':'六盘水市','盘州':'六盘水市','兴安盟':'兴安','大兴安岭地区':'大兴安岭'})
    city=city.str.replace(r'市$','',regex=True)
    counts=city.value_counts().rename_axis('City').reset_index(name='Subject_Count')
    assert counts.Subject_Count.sum()==len(culture)
    counts['Log_Subject_Count']=np.log(counts.Subject_Count)
    return counts

def generate(inputs,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    culture=pd.read_excel(inputs['culture_individual'])
    for field in ['province','city']:
        counts=city_counts(culture) if field=='city' else culture[field].value_counts().rename_axis(field.title()).reset_index(name='Subject_Count')
        counts['Log_Subject_Count']=np.log(counts.Subject_Count)
        counts.to_csv(output/f'Culture_{field}_counts.csv',index=False,encoding='utf-8-sig')
    residence=extract_primary_residence(inputs['residence_history'])
    cul=pd.read_excel(inputs['culture_grid'])
    c=cul.loc[cul.longitude.between(73.5,135,inclusive='neither')&cul.latitude.between(4,53,inclusive='neither')]
    im=residence.loc[residence.Primary_Longitude_WGS84.between(73.5,135,inclusive='neither')&residence.Primary_Latitude_WGS84.between(4,53,inclusive='neither')]
    cx,cy=lon_lat_2dikaer(c.longitude,c.latitude);ix,iy=lon_lat_2dikaer(im.Primary_Longitude_WGS84,im.Primary_Latitude_WGS84)
    wn_x,wn_y=lon_lat_2dikaer(80,60);en_x,_=lon_lat_2dikaer(140,60);_,ws_y=lon_lat_2dikaer(80,5)
    for step in [50,70,75,100,125,150,200,250,300]:
        xt=np.arange(wn_x,en_x+step,step);yt=np.arange(wn_y,ws_y+step,step)
        def counts(x,y):
            return pd.DataFrame(dict(x=get_matlab_grid_idx(xt,x),y=get_matlab_grid_idx(yt,y))).groupby(['x','y']).size()
        cc,ii=counts(cx,cy),counts(ix,iy);rows=[]
        for gx,gy in sorted(set(cc.index)|set(ii.index)):
            nc,ni=int(cc.get((gx,gy),0)),int(ii.get((gx,gy),0));x=xt[gx-1]+step/2;y=yt[gy-1]+step/2;lo,la=xy_to_coor(x,y)
            code=3 if nc and ni else (1 if nc else 2)
            rows.append([f'G{step}k_X{gx:03d}_Y{gy:03d}',gx+gy*1000,gx,gy,round(x,2),round(y,2),round(float(lo),6),round(float(la),6),nc,ni,nc+ni,code])
        pd.DataFrame(rows,columns=['Grid_ID','GroupID','xgrid','ygrid','Center_X_km','Center_Y_km','Center_Longitude_WGS84','Center_Latitude_WGS84','N_Cultural_Subjects','N_Imaging_Subjects','Total_Subjects','Category_Code']).to_csv(output/f'grid_{step}km.csv',index=False)
    print('Saved individual-derived counts, residence and grid tables',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    generate(json.loads(a.inputs.read_text(encoding='utf-8')),a.output)
