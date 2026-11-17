"""Create full-cohort country/site keys from the cortical input identifiers."""
import argparse
from pathlib import Path
import h5py,numpy as np,pandas as pd
from scipy.io import loadmat


def build(nations,appended,output):
    with h5py.File(nations) as h:
        ids=[''.join(chr(int(c)) for c in h[r][()].ravel()) for r in h['nationsdata_all/lh_500_aparc_thickness'][()].ravel()]
    ids=[s.replace('Turkey-1-','Mexico-3-').replace('Spain-2-','Belgium-3-') for s in ids if not s.startswith('Israel-') and not (s.startswith('Germany-3-') and 'acq-MP2R' in s)]
    us=loadmat(appended,simplify_cells=True)
    for j,state in enumerate(['Louisiana','Missouri','Wisconsin'],1):
        ids.extend(f'USA-{j}-{us["subid"][i]}' for i in np.flatnonzero(us['state']==state))
    aliases={'UK':'United Kingdom','USA':'United States'}
    table=pd.DataFrame(dict(row_key=np.arange(1,len(ids)+1),participant=ids,unit=[aliases.get(s.split('-')[0],s.split('-')[0]) for s in ids],batch=[s.split('-')[1] for s in ids]))
    output.parent.mkdir(parents=True,exist_ok=True);table.to_csv(output,index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['nations','appended','output']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();build(a.nations,a.appended,a.output)
