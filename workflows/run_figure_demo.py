"""Create synthetic inputs for the figure assembly examples."""
from pathlib import Path
import argparse
import numpy as np
from assemble import assemble,FIGURES
from panel_data import save
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'outputs/figure_demo');a=p.parse_args()
rng=np.random.default_rng(19);R=32;G=10;N=100;K=6
d=dict(effect_maps=rng.normal(size=(1,R)),gradients=rng.normal(size=(R,K)),null_maps=rng.normal(size=(24,1,R)),spin_indices=np.array([rng.permutation(R) for _ in range(24)]),cluster_influence=rng.normal(size=(1,G,R))*.03,
 province=np.repeat(np.arange(G),10),center=np.repeat(np.arange(G),10),coordinates_a_km=rng.normal(size=(N,2))*100,coordinates_b_km=rng.normal(size=(N,2))*100,cell_size_km=100,
 cortex=rng.normal(size=(N,R)),culture=np.repeat(rng.normal(size=(G,1)),10,axis=0),correlation_method=np.array('pearson'),
 measurements=rng.normal(size=(N,3)),measurement_labels=np.array(['A','B','C']),
 task_effects=rng.normal(size=(7,R)),best_order=np.arange(7),task_labels=np.array(list('ABCDEFG')),example_roi=np.array([0]),
 components=np.arange(5),statistics=rng.normal(size=(5,2)),series_labels=np.array(['A','B']),
 correlation_draws=rng.uniform(.5,1,(24,5,2)),proportions=np.linspace(.2,1,5),full_map=rng.normal(size=R),fold_maps=rng.normal(size=(G,R)),unit_labels=np.array([str(i) for i in range(G)]),
 estimate=np.arange(4)*.1,low=np.arange(4)*.1-.2,high=np.arange(4)*.1+.2,roi_labels=np.array(list('ABCD')),
 cultural_dimensions=rng.normal(size=(G,6)),dimension_labels=np.array(list('ABCDEF')),china_correlations=rng.normal(size=(3,2)),world_correlations=rng.normal(size=(3,2)),
 word_scores=rng.normal(size=(G,3)),validity_score=rng.normal(size=G),score_labels=np.array(list('ABC')),
 participant1=2.5+rng.normal(size=(G,R))*.1,participant2=2.5+rng.normal(size=(G,R))*.1,center_labels=np.array([str(i) for i in range(G)]))
for g in ['HH','LL','HL','LH','MHH','MLL']:d[g]=rng.normal(size=(20,2))

d.update(duration=rng.integers(0,10,N),birth_province=np.array(['P'+str(i%6) for i in range(N)]),destination_province=np.array(['P'+str(i%5) for i in range(N)]))
for key in ['SAaxis','mode2_display','xu_evo','reardon_dev','intervar','cogPC1','myelin','gene_PC1','fcgradient1']:d[key]=rng.normal(size=R)
d.update(china_dimension_effects=rng.normal(size=(R,6)),nations_dimension_effects=rng.normal(size=(R,6)),china_dimension_nulls=rng.normal(size=(24,R,6)),nations_dimension_nulls=rng.normal(size=(24,R,6)),gradient_labels=np.array(['Gradient'+str(i+1) for i in range(6)]))
for figure in FIGURES:
    inputs=d
    if figure=='Figure5':
        inputs=dict(d,webb_map_draws=d['effect_maps'][None]+rng.normal(size=(24,1,R))*.05)
    source,panels=assemble(figure,inputs,repetitions=24,seed=42)
    print(save(figure,a.output,source,panels))
