from pathlib import Path
import argparse,sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analyses'))
from word_scores import T_WORDS,L_WORDS
from run_specialized import run

DICTIONARY=np.array([['constrain','约束','allow','允许'],['uniformity','一致','freedom','自由'],['enforce','执行','create','创造'],['proscribe','禁止','variability','变化'],['abide','遵从','autonomy','自主'],['circumscribe','限制','openness','开放'],['uphold','维护','flexibility','灵活'],['forbid','不准','adapt','适应'],['strict','严格','encourage','鼓励'],['force','强制','choose','选择']],object)

def bridge(frame,output):
    words=T_WORDS+L_WORDS
    d=dict(word_counts=frame[words].to_numpy(float),scores=frame[['pc1_Diff_clean','pc1_T_clean','pc1_L_clean']].to_numpy(float),validity_score=frame.z_tightness_control.to_numpy(float),province_names=frame.names.to_numpy(str),word_labels=np.array(words,object),dictionary=DICTIONARY)
    return run('word_layout',d,output)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();bridge(pd.read_csv(a.input),a.output)
