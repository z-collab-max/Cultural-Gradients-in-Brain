from pathlib import Path
import argparse
from word_counts import generate as clean
from word_scores import generate as score

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--raw-input',type=Path,required=True)
    parser.add_argument('--province-input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    words=clean(args.raw_input,args.province_input,args.output/'cleaned')
    frame=score(words,args.province_input,args.output/'scores')
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
    from bridge_word_scores import bridge
    bridge(frame,args.output/'figure_source')
