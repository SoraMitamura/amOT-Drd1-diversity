"""Summarize all six dependent-correlation contrasts in one figure."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
from display_helpers import ranking_text

def plot_results(out):
    out=Path(out);d=pd.read_csv(out/'test_results.csv');rank=pd.read_csv(out/'similarity_ranking.csv')
    meta=json.loads((out/'manifest.json').read_text());title=meta['category'].removeprefix('GOBP_').replace('_',' ').capitalize()
    plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':13,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(14,9),facecolor='white')
    fig.text(.06,.96,'amOT · R | '+title,fontsize=20,weight='bold')
    fig.text(.06,.92,'Marker-excluded: '+(', '.join(meta['excluded_markers']) or 'none present'),fontsize=15,color='#a44335')
    fig.text(.06,.88,'Spearman ρ:  '+ranking_text(rank,d),fontsize=15)
    ax=fig.add_axes([.36,.32,.30,.44])
    y=np.arange(len(d))[::-1]
    for k,row in d.iterrows():
        ax.errorbar(row.delta_A_minus_B,y[k],xerr=[[row.delta_A_minus_B-row.ci95_low],[row.ci95_high-row.delta_A_minus_B]],fmt='o',color='#246a96',capsize=4,markersize=7,lw=2)
    ax.axvline(0,color='#777777',ls='--',lw=1)
    ax.set_yticks(y);ax.set_yticklabels([f'{r.A} − {r.B}' for r in d.itertuples()])
    ax.tick_params(axis='both',labelsize=17);ax.set_ylim(-.6,5.6);ax.set_xlabel('Δρ = ρ(tD1, A) − ρ(tD1, B)\nPoint: observed difference; bar: 95% percentile CI',fontsize=15)
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0,pad=14)
    for k,row in d.iterrows():
        yp=.32+.44*(y[k]+.6)/6.2
        fig.text(.71,yp,f'{row.p_unadjusted:.4g}',va='center',fontsize=17)
        fig.text(.82,yp,(row.closer_to_tD1 if row.p_unadjusted<.05 else 'n.s.'),va='center',fontsize=17)
    fig.text(.71,.82,'Raw p',weight='bold');fig.text(.82,.82,'Closer to tD1',weight='bold',fontsize=15)
    fig.text(.06,.215,'≈ (n.s.): no significant difference (raw p ≥ 0.05).',fontsize=14)
    fig.text(.06,.175,'Positive Δρ: A is more similar to tD1 than B.',fontsize=14,weight='bold')
    B=int(d.resamples.iloc[0]);floor=1/(B+1)
    fig.text(.06,.125,f'Two-sided centered cell bootstrap; {B:,} resamples; {meta["genes"]} genes; no section stratification.',fontsize=15)
    fig.text(.06,.085,f'Raw p-values shown. Monte Carlo resolution floor: {floor:.4g}. Holm-adjusted values are included in CSV.',fontsize=15)
    fig.text(.06,.045,'Conditional on one specimen’s fixed imputed data; cells are not independent biological replicates.',fontsize=13,color='#505b65')
    for ext in ['png','svg','pdf']:fig.savefig(out/f'all_six_correlation_tests.{ext}',dpi=300)
    plt.close(fig)
