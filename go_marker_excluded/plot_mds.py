"""Reuse the prior Spearman chord-distance metric MDS on the selected GO set."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist,squareform
from sklearn.manifold import MDS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from data_io import GROUPS

def plot_mds(out):
    out=Path(out);names=[n for n,_ in GROUPS]
    meta=json.loads((out/'manifest.json').read_text());title=meta['category'].removeprefix('GOBP_').replace('_',' ').capitalize()
    df=pd.read_csv(out/'gene_means_and_ranks.csv')
    x=df[[f'mean_{n}' for n in names]].to_numpy().T
    corr=spearmanr(x,axis=1).statistic
    dist=np.sqrt(np.maximum(0,2*(1-corr)));np.fill_diagonal(dist,0)
    assert np.isfinite(dist).all() and np.allclose(dist,dist.T)
    xy=MDS(n_components=2,metric=True,dissimilarity='precomputed',random_state=0,n_init=32,max_iter=3000,eps=1e-9).fit_transform(dist)
    projected=squareform(pdist(xy));upper=np.triu_indices(5,1)
    stress=float(np.sqrt(np.sum((projected[upper]-dist[upper])**2)/np.sum(dist[upper]**2)))
    for a,name in [(corr,'spearman_correlations'),(dist,'correlation_distances')]:pd.DataFrame(a,index=names,columns=names).to_csv(out/f'{name}.csv')
    pd.DataFrame(xy,index=names,columns=['MDS1','MDS2']).to_csv(out/'MDS_coordinates.csv')
    pd.DataFrame([dict(A=names[i],B=names[j],rho=corr[i,j],original_distance=dist[i,j],distance_2D=projected[i,j]) for i,j in zip(*upper)]).to_csv(out/'distance_reconstruction.csv',index=False)
    plt.rcParams.update({'font.family':'Arial','font.size':14,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(6.2,6.2))
    colors=['#D94732','#E69F00','#0072B2','#009E73','#8454B3']
    for pt,n,c in zip(xy,names,colors):
        ax.scatter(*pt,s=100,c=c,edgecolors='white',zorder=3)
    ax.set_aspect('equal',adjustable='box');ax.margins(.3)
    center=(xy.max(axis=0)+xy.min(axis=0))/2
    half=float(np.ptp(xy,axis=0).max())*.78
    ax.set_xlim(center[0]-half,center[0]+half)
    ax.set_ylim(center[1]-half,center[1]+half)
    ax.set_xlabel('MDS 1');ax.set_ylabel('MDS 2');ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
    import textwrap
    ax.set_title('amOT · R | Marker-excluded\n'+textwrap.fill(title,38),fontsize=14,pad=16)
    fig.subplots_adjust(left=.16,right=.96,top=.80,bottom=.32)
    removed=', '.join(meta['excluded_markers']) or 'none present'
    fig.text(.5,.18,'Excluded: '+removed,ha='center',fontsize=11,color='#a44335')
    fig.text(.5,.12,f'{len(df)} genes | Distance = √[2(1 − Spearman ρ)]',ha='center',fontsize=11)
    fig.text(.5,.055,f'Metric MDS | normalized distance stress = {stress:.3f}',ha='center',fontsize=11)
    from display_helpers import place_labels
    place_labels(fig,ax,xy,names,colors)
    for ext in ['png','svg','pdf']:fig.savefig(out/f'GO_guidance_MDS.{ext}',dpi=300)
    plt.close(fig)
    manifest=json.loads((out/'manifest.json').read_text())
    manifest.update(MDS_stress=stress,MDS_seed=0,MDS_n_init=32)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('MDS stress:',stress)

if __name__=='__main__':plot_mds(Path(__file__).parent/'results')
