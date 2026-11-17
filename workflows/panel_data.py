from pathlib import Path
import json,re
import numpy as np
from scipy.io import savemat


def panel(figure,name,kind,**values):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name):raise ValueError('Invalid panel name')
    return dict(folder=figure,name=name,kind=kind,**values)


def labels(prefix,n):return np.array([f'{prefix}{j+1}' for j in range(n)],dtype=object)


def save(figure,output,source,panels):
    output=Path(output)/figure;output.mkdir(parents=True,exist_ok=True)
    safe={k:np.asarray(v).astype(str) if np.asarray(v).dtype.kind=='O' else v for k,v in source.items()}
    np.savez_compressed(output/'source_data.npz',**safe)
    savemat(output/'source_data.mat',source,do_compression=True)
    cells=np.empty((1,len(panels)),object)
    for j,p in enumerate(panels):cells[0,j]=p
    savemat(output/'panels.mat',dict(panels=cells),do_compression=True)
    schema=[dict(name=p['name'],kind=p['kind'],fields=sorted(p)) for p in panels]
    (output/'panels.json').write_text(json.dumps(schema,indent=2),encoding='utf-8')
    return output
