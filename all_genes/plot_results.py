"""Summarize all six dependent-correlation contrasts in one figure."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_results(out):
    out=Path(out);d=pd.read_csv(out/'test_results.csv');rank=pd.read_csv(out/'similarity_ranking.csv')
    plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':13,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(12,8),facecolor='white')
    fig.text(.06,.94,'amOT · R | Which population is more similar to tD1?',fontsize=22,weight='bold')
    fig.text(.06,.88,'Spearman ρ:  '+'   >   '.join(f'{r.comparator} ({r.rho_with_tD1:.4f})' for r in rank.itertuples()),fontsize=13)
    ax=fig.add_axes([.34,.31,.38,.48])
    y=np.arange(len(d))[::-1]
    for k,row in d.iterrows():
        ax.errorbar(row.delta_A_minus_B,y[k],xerr=[[row.delta_A_minus_B-row.ci95_low],[row.ci95_high-row.delta_A_minus_B]],fmt='o',color='#246a96',capsize=4,markersize=7,lw=2)
    ax.axvline(0,color='#777777',ls='--',lw=1)
    ax.set_yticks(y);ax.set_yticklabels([f'{r.A} − {r.B}' for r in d.itertuples()])
    ax.set_ylim(-.6,5.6);ax.set_xlabel('Δρ = ρ(tD1, A) − ρ(tD1, B)\nPoint: observed difference; bar: 95% percentile CI',fontsize=12)
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0,pad=14)
    for k,row in d.iterrows():
        yp=.31+.48*(y[k]+.6)/6.2
        fig.text(.77,yp,f'{row.p_unadjusted:.4g}',va='center',fontsize=13)
        fig.text(.86,yp,row.closer_to_tD1,va='center',fontsize=12)
    fig.text(.77,.82,'Raw p',weight='bold');fig.text(.86,.82,'Closer to tD1',weight='bold',fontsize=12)
    fig.text(.06,.20,'Positive Δρ: A is more similar to tD1 than B.',fontsize=14,weight='bold')
    B=int(d.resamples.iloc[0]);floor=1/(B+1)
    fig.text(.06,.14,f'Two-sided centered cell bootstrap; {B:,} resamples; 8,460 genes; no section stratification.',fontsize=12)
    fig.text(.06,.10,f'Raw p-values shown. Monte Carlo resolution floor: {floor:.4g}. Holm-adjusted values are included in CSV.',fontsize=11)
    fig.text(.06,.06,'Conditional on one specimen’s fixed imputed data; cells are not independent biological replicates.',fontsize=11,color='#505b65')
    for ext in ['png','svg','pdf']:fig.savefig(out/f'all_six_correlation_tests.{ext}',dpi=300)
    plt.close(fig)
