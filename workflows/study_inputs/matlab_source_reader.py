"""Read numeric arrays, cells and structs from MATLAB source files."""
import h5py
import numpy as np
from scipy.io import loadmat

def read_node(h, obj):
    if isinstance(obj, h5py.Group):
        return {k: read_node(h, v) for k, v in obj.items()}
    if obj.attrs.get('MATLAB_empty', 0):
        return np.array([])
    a = obj[()]
    if h5py.check_dtype(ref=a.dtype):
        values = [read_node(h, h[r]) for r in a.T.ravel()]
        return values[0] if len(values) == 1 else values
    if obj.attrs.get('MATLAB_class') == b'char':
        return ''.join(chr(int(n)) for n in a.ravel())
    return a.T.squeeze()

def read_mat(path, key):
    if h5py.is_hdf5(path):
        with h5py.File(path) as h:
            return read_node(h, h[key])
    return loadmat(path, simplify_cells=True)[key]
