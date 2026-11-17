"""Prepare model tables from raw macro variables at their geographic level."""
import argparse,json
from pathlib import Path
import pandas as pd
import numpy as np
from macro_pca import fit_pca
from input_fields import CHINA_CITY,CHINA_PROVINCE,NATIONS_BLOCKS,USA_MACROS


def prepare(table, blocks):
    out=table.copy();metadata={}
    for prefix,keys,fields in blocks:
        if any(str(v).startswith(prefix+'_pc') for v in table):raise ValueError('Use raw inputs without existing '+prefix+' PCs')
        if not set(keys+fields)<=set(table):raise ValueError('Missing PCA fields: '+prefix)
        if table[keys+fields].isna().any().any():raise ValueError('Missing geographic keys or macro values')
        grouped=table.groupby(keys,sort=True,dropna=False)
        if (grouped[fields].nunique()>1).any().any():
            raise ValueError('Macro variables must be constant within each '+prefix+' unit')
        units=grouped[fields].first().reset_index()
        scores,meta=fit_pca(units[fields].to_numpy(float),prefix)
        meta.update(fields=fields,keys=keys,unit_keys=units[keys].astype(str).values.tolist())
        pc=units[keys].copy()
        for j,name in enumerate(meta['names']):pc[name]=scores[:,j]
        out=out.drop(columns=fields).merge(pc,on=keys,how='left',validate='many_to_one',sort=False)
        metadata[prefix]=meta
    return out,metadata


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort',choices=['China','Nations','USA'],required=True)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.input.resolve()==a.output.resolve():raise ValueError('Use a separate output file')
    keys={'China':['participant_id','province_birth','city_birth'],'Nations':['unit'],'USA':['unit']}[a.cohort]
    table=pd.read_csv(a.input,dtype={k:str for k in keys})
    if a.cohort=='China':
        if table.participant_id.isna().any() or table.participant_id.duplicated().any():raise ValueError('Unique participant IDs required')
        blocks=[('city_macro',['province_birth','city_birth'],CHINA_CITY),('province_macro',['province_birth'],CHINA_PROVINCE)]
    else:
        if table.unit.isna().any() or table.unit.duplicated().any():raise ValueError('Supply one row per geographic unit')
        blocks=[('nation_macro' if a.cohort=='Nations' else 'state_macro',['unit'],sum(NATIONS_BLOCKS.values(),[]) if a.cohort=='Nations' else USA_MACROS)]
    out,meta=prepare(table,blocks);a.output.parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.output,index=False)
    a.output.with_suffix('.pca.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')

if __name__=='__main__':main()
