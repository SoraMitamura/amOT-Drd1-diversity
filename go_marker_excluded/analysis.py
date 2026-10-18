"""Six dependent Spearman correlation contrasts, sharing tD1 in every replicate."""
from pathlib import Path
from itertools import combinations
import argparse,json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
from threadpoolctl import threadpool_limits
from data_io import GROUPS,load_profile,gene_key

def corr_rows(a,b):
    """Pearson correlation of paired rank rows (Spearman); genes are not resampled."""
    a=a-a.mean(axis=1,keepdims=True);b=b-b.mean(axis=1,keepdims=True)
    return np.sum(a*b,axis=1)/np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b,axis=1))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,default=Path('results'))
    parser.add_argument('--resamples',type=int,default=9999)
    parser.add_argument('--seed',type=int,default=20261017)
    parser.add_argument('--selection',type=Path,required=True)
    args=parser.parse_args()
    if args.resamples<100:parser.error('At least 100 resamples required')
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    matrices=[];headers=[];sources=[]
    for name,prefix in GROUPS:
        path=args.input_dir/f'{prefix}_OT_R_amOT_v3.csv'
        m,genes,n=load_profile(path);matrices.append(m);headers.append(genes)
        sources.append(dict(group=name,file=path.name,n_cells=n,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        print(f'{name}: {n} cells',flush=True)
    genes=headers[1];keys=[gene_key(g) for g in genes]
    if any(set(m.columns)!=set(keys) for m in matrices):raise ValueError('Gene sets differ')
    selection=json.loads(args.selection.read_text(encoding='utf-8'))
    genes=[r['gene'] for r in selection['records'] if r['in_expression_panel']]
    keys=[gene_key(g) for g in genes]
    assert len(keys)>=10 and len(set(keys))==len(keys)
    matrices=[m[keys].to_numpy() for m in matrices]
    names=[g[0] for g in GROUPS];means=np.stack([m.mean(axis=0) for m in matrices])
    ranks=rankdata(means,axis=1)
    observed=np.array([spearmanr(means[0],means[i]).statistic for i in range(1,5)])
    assert np.allclose(observed,[corr_rows(ranks[:1],ranks[i:i+1])[0] for i in range(1,5)])
    # A synthetic invariant checks the shared-reference correlation calculation.
    a=np.array([[1.,2.,4.,3.],[4.,3.,1.,2.]])
    assert np.allclose(corr_rows(a,a),1) and np.allclose(corr_rows(a,-a),-1)
    profile=pd.DataFrame({'gene':genes})
    for i,name in enumerate(names):profile[f'mean_{name}']=means[i];profile[f'rank_{name}']=ranks[i]
    profile.to_csv(out/'gene_means_and_ranks.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(sources).to_csv(out/'input_provenance.csv',index=False,encoding='utf-8-sig')
    rng=np.random.default_rng(args.seed);B=args.resamples;boots=np.empty((B,4));batch=64
    with threadpool_limits(limits=4):
        for start in range(0,B,batch):
            k=min(batch,B-start)
            profiles=[rng.multinomial(len(x),np.full(len(x),1/len(x)),size=k).astype(float)@x/len(x) for x in matrices]
            ranked=[rankdata(p,axis=1) for p in profiles]
            for j in range(4):boots[start:start+k,j]=corr_rows(ranked[0],ranked[j+1])
            if start%1024==0:print(f'Bootstrap {start+k}/{B}',flush=True)
    bootdf=pd.DataFrame(boots,columns=[f'rho_tD1_{n}' for n in names[1:]])
    bootdf.insert(0,'replicate',np.arange(1,B+1))
    rows=[]
    for i,j in combinations(range(4),2):
        delta=observed[i]-observed[j];ds=boots[:,i]-boots[:,j];null=ds-delta
        tail=np.abs(null)>=abs(delta);low,high=np.quantile(ds,[.025,.975]);p=(1+tail.sum())/(B+1)
        tag=f'{names[i+1]}_minus_{names[j+1]}'
        bootdf[f'delta_{tag}']=ds;bootdf[f'centered_null_{tag}']=null;bootdf[f'tail_exceedance_{tag}']=tail
        rows.append(dict(reference='tD1',A=names[i+1],B=names[j+1],rho_tD1_A=observed[i],rho_tD1_B=observed[j],delta_A_minus_B=delta,ci95_low=low,ci95_high=high,p_unadjusted=float(p),null_tail_count=int(tail.sum()),resamples=B,closer_to_tD1=names[i+1] if delta>0 else names[j+1]))
    results=pd.DataFrame(rows)
    # Holm adjustment supplied transparently as an additional six-comparison check.
    order=np.argsort(results.p_unadjusted.to_numpy());adjusted=np.empty(6)
    adjusted[order]=np.minimum(1,np.maximum.accumulate(results.p_unadjusted.to_numpy()[order]*(6-np.arange(6))))
    results['p_holm_six_comparisons']=adjusted
    results.to_csv(out/'test_results.csv',index=False,encoding='utf-8-sig')
    bootdf.to_csv(out/'bootstrap_calculations.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame({'comparator':names[1:],'rho_with_tD1':observed}).sort_values('rho_with_tD1',ascending=False).to_csv(out/'similarity_ranking.csv',index=False,encoding='utf-8-sig')
    steps=[('1','Within each group, average supplied log2 values across cells; all sections pooled.'),('2','Rank the 196 available GO-category genes within each mean profile; average ranks for ties.'),('3','Compute four Spearman correlations with tD1; contrast delta=rho(tD1,A)-rho(tD1,B).'),('4','Independently resample whole cells with replacement within each group, retaining its original size; reuse the SAME tD1 sample for all four correlations in each replicate.'),('5','Recompute means, ranks, correlations and all six contrasts in each replicate.'),('6','95% percentile CI = quantile(delta_star,[0.025,0.975]).'),('7','Approximate null = delta_star - observed_delta. Two-sided p=(1+count(abs(null)>=abs(observed_delta)))/(B+1).'),('8','Report unadjusted p primarily; Holm-adjusted p also provided for the family of six comparisons.'),('9','Cell-exchangeability assumption; conditional on fixed imputed values. Spatial dependence and imputation uncertainty are not propagated. Not animal-level inference.')]
    pd.DataFrame(steps,columns=['step','description']).to_csv(out/'calculation_steps.csv',index=False,encoding='utf-8-sig')
    (out/'manifest.json').write_text(json.dumps(dict(sources=sources,genes=len(genes),seed=args.seed,resamples=B,batch_size=batch,anatomy='right amOT v3, all sections pooled',n_specimens=1,notes='Exploratory centered cell bootstrap, not a label permutation test or animal-level test. Entire five-group bootstrap freshly rerun; Monte Carlo values can differ slightly from prior three-group runs.'),indent=2),encoding='utf-8')
    from plot_results import plot_results
    manifest=json.loads((out/'manifest.json').read_text())
    manifest.update(category=selection['category'],GO_total=selection['total'],available_genes=len(genes),missing_genes=selection['missing'],excluded_markers=selection['excluded_markers'],exclusion_panel=selection['exclusion_panel'],analysis='Marker-excluded sensitivity analysis; original GO definition retained in source, analysis subset modified explicitly.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    steps_path=out/'calculation_steps.csv'
    steps_path.write_text(steps_path.read_text(encoding='utf-8-sig').replace('196 available GO-category genes',f'{len(genes)} available GO-category genes'),encoding='utf-8-sig')
    plot_results(out)
    from plot_mds import plot_mds
    plot_mds(out)
    print(results.to_string(index=False),flush=True)

if __name__=='__main__':main()
