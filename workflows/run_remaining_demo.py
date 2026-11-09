"""Synthetic checks for the specialized source-data adapters."""
from pathlib import Path
import argparse,json,sys
import numpy as np
from run_specialized import run,FIGURES
from bridge_word_scores import DICTIONARY,T_WORDS,L_WORDS

def main(output):
    output=Path(output);rng=np.random.default_rng(20261103);N=72;R=12;G=6;B=24
    province=np.repeat(np.arange(G),N//G);grad=rng.normal(size=(R,6));full=rng.normal(size=(2,R))
    maps=dict(effect_maps=full,gradients=grad,null_maps=rng.normal(size=(B,2,R)),spin_indices=np.array([rng.permutation(R) for _ in range(B)]),left_hemisphere=np.arange(R)<R//2,model_labels=np.array(['Model_A','Model_B']))
    base=dict(x=rng.normal(size=(N,2)),y=rng.normal(size=(N,2)),province=province)
    inputs={'province_summary':base,'ecological_comparison':base,
       'demographics':dict(continuous=rng.normal(size=(N,2)),continuous_labels=np.array(['Age','Education']),categorical=rng.integers(0,2,N),group=province),
       'gradients':{key:rng.normal(size=R) for key in ['SAaxis','mode2_display','xu_evo','reardon_dev','intervar','cogPC1','myelin','gene_PC1','fcgradient1']},
       'commonality':dict(mover_maps=rng.normal(size=(2,R)),stayer_maps=rng.normal(size=(2,R)),deculturation_maps=rng.normal(size=(2,R)),acculturation_maps=rng.normal(size=(2,R)),relocation_maps=rng.normal(size=(2,R))),
       'hierarchy':dict(task_effects=rng.normal(size=(3,R,7)),branch_labels=np.array(['Behavior','PC1','PC2']),task_labels=np.array(list('ABCDEFG')),example_roi=np.array([0,1]),spin_indices=maps['spin_indices'],left_hemisphere=maps['left_hemisphere']),
       'rice_culture':dict(group_values=rng.normal(size=(30,4)),labels=np.array(['RiceShare','RiceSuitability','PC1','PC2'])),
       'effect_pca':dict(dimension_effects=rng.normal(size=(R,6)),culture_effects=full.T,reference_component=np.array([0,1]),dimension_labels=np.array(['D'+str(i) for i in range(6)]),spin_indices=maps['spin_indices'],left_hemisphere=maps['left_hemisphere']),
       'resampling':dict(full_maps=full,refit_maps=full[:,None,None,:]+rng.normal(size=(2,B,4,R))*.2,proportions=np.array([.25,.5,.75,1]),series_labels=np.array(['Group','Individual'])),
       'traveling_unshared':dict(participant1=2.5+rng.normal(size=(11,32))*.1,participant2=2.5+rng.normal(size=(11,32))*.1,centers1=np.array(['S'+str(i) for i in range(11)]),centers2=np.array(['S'+str(i) for i in range(2,13)])),
       'network':{g:rng.normal(size=(20,2)) for g in ['HH','LL','HL','LH','MHH','MLL']},
       'reference_histograms':dict(effect_maps=full,null_maps=maps['null_maps'],reference=grad[:,0]),
       'dimension_omission':dict(full_maps=full,fold_maps=full[None,:,:]+rng.normal(size=(6,2,R))*.1,dimension_labels=np.array(['D'+str(i) for i in range(6)]),spin_indices=maps['spin_indices']),
       'thresholds':dict(endpoint_map=full[0],threshold_maps=full[0]+rng.normal(size=(4,R))*.2,threshold_labels=np.array(['90%','95%','99%','endpoint'])),
       'overlap_curves':dict(scales=np.array([50,100,200]),overlap_counts=np.array([[10,8,7],[9,5,4],[8,2,1]])),
       'psychometrics':dict(dimension_scores=rng.normal(size=(N,2)),items=rng.integers(1,6,(N,8)),item_mask=np.array([[1]*4+[0]*4,[0]*4+[1]*4],bool),response_max=np.array([5,5]),compute_alpha=np.array([False,True]),anova_group=province,icc_group=province,dimension_labels=np.array(['Dimension','TL'])),
       'brain_icc':dict(adjusted_brain=rng.normal(size=(3,N,R)),province=province,measure_labels=np.array(['Area','Thickness','Volume']))}
    for name in ['combat_versions','area_volume','rice_brain']:inputs[name]=maps.copy()
    inputs['international']=dict(maps,webb_map_draws=maps['effect_maps'][None]+rng.normal(size=(B,2,R))*.02)
    tight=rng.normal(size=31);loose=rng.normal(size=31)
    inputs['word_layout']=dict(word_counts=rng.integers(20,500,(31,20)),scores=np.c_[tight-loose,tight,loose],validity_score=rng.normal(size=31)*.05,province_names=np.array(['P'+str(i) for i in range(31)]),word_labels=np.array(T_WORDS+L_WORDS),dictionary=DICTIONARY.astype(str))
    assert set(inputs)==set(FIGURES)
    schema={};roots=[]
    for name,d in inputs.items():
        source=output/'inputs';source.mkdir(parents=True,exist_ok=True);np.savez_compressed(source/(name+'.npz'),**d)
        roots.append(str(run(name,d,output,B,42).relative_to(output)))
        schema[name]={k:dict(shape=list(np.asarray(v).shape),dtype=str(np.asarray(v).dtype)) for k,v in d.items()}
    (output/'run_index.json').write_text(json.dumps(roots,indent=2),encoding='utf-8')
    (output/'input_schema.json').write_text(json.dumps(schema,indent=2),encoding='utf-8')
    print('SPECIALIZED_DEMO_PASS',len(roots))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'outputs/specialized');a=p.parse_args();main(a.output)
