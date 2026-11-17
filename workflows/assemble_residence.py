"""Build SI4 source data from aligned annual histories and participant metadata."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from assemble import assemble
from panel_data import save
from scipy.io import savemat


def build(history, participants, output):
    table = pd.read_csv(participants, dtype={'id': 'int64'})
    with np.load(history, allow_pickle=False) as records:
        if table.id.duplicated().any() or not np.array_equal(table.id, records['ids']):
            raise ValueError('History IDs must match participant rows in order')
        selected = table.valid.astype(str).str.lower().isin(['1', 'true']) & table.lived.eq(2)
        routes = table.loc[selected, ['birth', 'destination']]
        if routes.isna().any().any():
            raise ValueError('Selected routes need both province labels')
        names = json.loads(Path(__file__).with_name('province_names.json').read_text(encoding='utf-8-sig'))
        routes = routes.replace(names)
        inputs = dict(duration=records['duration'][selected], birth_province=routes.birth.to_numpy(str),
                      destination_province=routes.destination.to_numpy(str))
    source, panels = assemble('SupplementaryFigure4', inputs)
    folder = save('SupplementaryFigure4', output, source, panels)
    # MATLAB expects one label per cell, including Unicode province names.
    labels = np.empty(len(source['route_labels']), dtype=object)
    for i, label in enumerate(source['route_labels']):
        labels[i] = np.frombuffer(str(label).encode('utf-16-le'), dtype='<u2')
    savemat(folder / 'source_data.mat', dict(source, route_labels=labels))
    jobs = [dict(kind='layout', function='reproduce_suppfig4', input='source_data.mat')]
    (folder / 'render_manifest.json').write_text(json.dumps(jobs), encoding='utf-8')
    return folder


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', type=Path, required=True)
    parser.add_argument('--participants', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(build(args.history, args.participants, args.output))
