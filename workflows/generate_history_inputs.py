"""Create residence duration and route inputs from annual residence records."""
from pathlib import Path
import argparse,csv,json,zipfile,xml.etree.ElementTree as ET
from collections import defaultdict
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--history',type=Path,required=True)
parser.add_argument('--participants',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
participants=list(csv.DictReader(args.participants.open(encoding='utf-8-sig')))
hist=defaultdict(list);ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(args.history) as book:
    shared=[]
    if 'xl/sharedStrings.xml' in book.namelist():
        shared=[''.join(x.itertext()) for x in ET.fromstring(book.read('xl/sharedStrings.xml')).findall('s:si',ns)]
    with book.open('xl/worksheets/sheet1.xml') as fh:
        for _,row in ET.iterparse(fh,events=['end']):
            if row.tag!='{'+ns['s']+'}row':continue
            if row.get('r')=='1':row.clear();continue
            cells={}
            for cell in row.findall('s:c',ns):
                col=''.join(c for c in cell.get('r') if c.isalpha());value=cell.find('s:v',ns)
                if cell.get('t')=='inlineStr':val=''.join(cell.find('s:is',ns).itertext())
                elif value is None:continue
                elif cell.get('t')=='s':val=shared[int(value.text)]
                else:val=value.text
                cells[col]=val
            if 'A' in cells and float(cells['A']).is_integer():
                hist[int(float(cells['B']))].append((int(float(cells['C'])),str(cells['F'])))
            row.clear()
ids=np.array([int(float(x['id'])) for x in participants]);count=np.full(len(ids),np.nan)
duration=count.copy();move_age=count.copy();last_duration=count.copy();birth=np.full(len(ids),'',dtype='<U50')
for j,ident in enumerate(ids):
    if ident not in hist:continue
    byyear=defaultdict(set)
    for yr,prov in sorted(hist[ident]):byyear[yr].add(prov)
    if any(len(x)>1 for x in byyear.values()):continue
    h=sorted((yr,next(iter(v))) for yr,v in byyear.items());birth[j]=h[0][1]
    trans=[h[k][0] for k in range(1,len(h)) if h[k][1]!=h[k-1][1]];count[j]=len(trans)
    if trans:duration[j]=h[-1][0]-trans[0];move_age[j]=trans[0]-h[0][0];last_duration[j]=h[-1][0]-trans[-1]
np.savez_compressed(args.output/'history_features.npz',ids=ids,count=count,duration=duration,
    move_age=move_age,last_duration=last_duration,birth=birth)
valid=[r for r in participants if str(r['valid']).lower() in ['1','true']]
with (args.output/'participant_clusters.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['original_row','birth_province','destination_province'])
    w.writerows((r['row'],r['birth'],r['destination']) for r in valid)
print('Saved history features for',len(ids),'participants and',len(valid),'eligible route records.')
