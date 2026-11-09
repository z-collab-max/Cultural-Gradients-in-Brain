import numpy as np
import h5py
RICE_MODELS = {
    "PercentPaddy": "talhelm_percent_paddy",
    "RiceSuitability": "talhelm_rice_suit",
}
PREDICTOR_NAMES = ("CulturePC1", "CulturePC2", "Rice")
MAP_FIELDS = ("SAaxis", "reardon_dev", "xu_evo", "intervar", "cogPC1", "mode2")
MAP_NAMES = ("SAaxis", "devexp", "evoexp", "intersubjvar", "cogpc1", "eigenmode2")

def numeric(group: h5py.Group, name: str) -> np.ndarray:
    return np.asarray(group[name], dtype=float).reshape(-1)


def matlab_cellstr(mat_file: h5py.File, dataset: h5py.Dataset) -> list[str]:
    values = []
    for ref in np.asarray(dataset).reshape(-1):
        chars = np.asarray(mat_file[ref]).reshape(-1)
        values.append("".join(chr(int(x)) for x in chars))
    return values


def zcols(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = x.std(axis=0, ddof=1)
    keep = np.isfinite(sd) & (sd > 1e-12)
    return (x[:, keep] - x[:, keep].mean(axis=0)) / sd[keep]


def zrows(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = x.std(axis=1, ddof=1, keepdims=True)
    return (x - x.mean(axis=1, keepdims=True)) / sd


def zvec(x: np.ndarray) -> np.ndarray:
    return (x - x.mean()) / x.std(ddof=1)


