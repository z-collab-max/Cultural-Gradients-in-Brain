from pathlib import Path
import numpy as np
from scipy.io import savemat


def save_source(output,stem,**arrays):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/(stem+'.npz'),**arrays)
    savemat(output/(stem+'.mat'),arrays)


def load_npz(path,required):
    with np.load(path,allow_pickle=False) as archive:
        missing=set(required)-set(archive.files)
        if missing:raise ValueError('Missing input fields: '+', '.join(sorted(missing)))
        return {k:archive[k] for k in required}
