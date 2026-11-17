"""Assemble SI2 from left/right parcel-thickness exports for two participants."""
import argparse
from pathlib import Path

import pandas as pd
from run_specialized import run


def read_hemisphere(path):
    raw = pd.read_csv(path, sep='\t')
    centers = raw.iloc[:, 0].str.extract(r'_(C\d+(?:_[A-Za-z]+)?)\.', expand=False)
    if centers.isna().any() or centers.duplicated().any():
        raise ValueError('Each scan needs a distinct center identifier')
    columns = [name for name in raw if name.endswith('_thickness') and 'unknown' not in name]
    if not columns:
        raise ValueError('No parcel thickness columns found')
    return raw[columns].set_axis(centers, axis=0)


def build(files, output):
    inputs = {}
    reference = None
    for i, (left, right) in enumerate(files, start=1):
        table = pd.concat([read_hemisphere(left), read_hemisphere(right)], axis=1, join='inner')
        if reference is None:
            reference = table.columns
        if not reference.equals(table.columns):
            raise ValueError('Parcel columns and order must match between participants')
        inputs[f'participant{i}'] = table.to_numpy(float)
        inputs[f'centers{i}'] = table.index.to_numpy(str)
    return run('traveling_unshared', inputs, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=Path, nargs=2, required=True, metavar=('LH', 'RH'))
    parser.add_argument('--second', type=Path, nargs=2, required=True, metavar=('LH', 'RH'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(build([args.first, args.second], args.output))
