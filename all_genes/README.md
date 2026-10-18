# All six tD1-reference correlation contrasts

Compare `rho(tD1, A)` versus `rho(tD1, B)` for every pair of comparators among **D2, Sema5a-D1, pDwarf Cell, ICj**. This is six contrasts, not ten pairwise tests between all populations. Analyses use **right amOT v3**, pooling all available sections, and all 8,460 genes.

## Run

Python 3.11 or 3.12 recommended:

```bash
python -m pip install -r requirements.txt
python analysis.py --input-dir data --output-dir results
```

Place five files directly in `data/`: `c61_OT_R_amOT_v3.csv`, `c62_OT_R_amOT_v3.csv`, `c63_OT_R_amOT_v3.csv`, `c264_OT_R_amOT_v3.csv`, `c262_OT_R_amOT_v3.csv`. These are the original amOT datasets, **not** the recently supplied section-specific coronal datasets. The display label pDwarf Cell denotes the supplied c264 population.

## Input checks

The loader requires nonempty files, `partition_domain=amOT`, `partition_hemisphere=R`, `z_ccf<5.7`, and unique nonmissing `partition_source_row` within each file. It excludes known metadata and all `partition_` columns and requires exactly 8,460 gene columns. Missing or infinite expression is rejected. MATLAB-style gene names are reconciled through a bijection (punctuation to underscore, numeric-leading symbols prefixed by x); D2 supplies canonical gene symbols. Gene sets must match exactly. Cell labels are not deduplicated because previous exports contain rounded IDs.

## Statistical procedure

1. Average the supplied log2 values over cells within each group, without additional normalization, transformation, or rounding.
2. Rank the 8,460 gene means within each group, using average ranks for ties. Spearman correlation is Pearson correlation of these ranks.
3. Define each observed contrast as `delta = rho(tD1,A) - rho(tD1,B)`. A positive contrast favors A.
4. Independently sample whole cell vectors with replacement within each of the five populations, retaining each population's original number of cells. This preserves within-cell relationships across genes. Reuse the **same sampled tD1 profile** in all four correlations and all six contrasts for that replicate. Recompute means, ranks, and correlations.
5. Repeat 9,999 times (seed 20261017, batches of 64, four BLAS threads). Sampling uses multinomial weights, equivalent to cell sampling with replacement. Changing group order or batch size changes the random stream.
6. Calculate a 95% percentile confidence interval from the bootstrap contrasts.
7. Approximate the null by centering bootstrap contrasts: `null_delta = delta_star - delta_observed`. Report the two-sided p-value `(1 + count(abs(null_delta) >= abs(delta_observed))) / (B + 1)`.

This is a **centered bootstrap approximation**, not an exact permutation test, not a gene-independent Fisher-z test, and not a test of MDS distances. The Monte Carlo resolution floor is 0.0001 for 9,999 resamples; zero tail exceedances do not imply a true p-value of zero or exactly 0.0001.

Unadjusted p-values are primary and shown in the PNG as requested. Holm-adjusted p-values for the six-test family are additionally included in CSV. The bootstrap is rerun jointly for all five populations; observed correlations match previous analyses, while Monte Carlo intervals may differ slightly from older three-group runs.

## Interpretation

These are exploratory comparisons conditional on a **single specimen's fixed imputed matrix**. The cell bootstrap assumes exchangeability within each population. Spatial dependence and uncertainty from imputation are not propagated, and cells are not independent animal replicates. Do not interpret small p-values as demonstrated animal-level population differences. No section stratification is performed.

## Output

- `all_six_correlation_tests.png` / SVG / PDF: all six contrasts with confidence intervals and raw p-values
- `test_results.csv`: correlations, contrasts, confidence intervals, raw and Holm p-values, tail counts and direction
- `similarity_ranking.csv`: four correlations with tD1, descending
- `gene_means_and_ranks.csv`: values underlying observed correlations
- `bootstrap_calculations.csv`: all four correlations and six contrasts, centered values and tail indicators for every replicate
- `calculation_steps.csv`: calculation summary
- `input_provenance.csv`: input hashes and cell counts
- `manifest.json`: seed, parameters, interpretation and specimen count

Dependencies and code are portable. Input data and results are excluded from Git. `data_io.py` adapts the validated loader from the prior MDS/scatter package to return matrices. `corr_rows` extends the prior correlation-contrast workflow to all five groups and checks identity/inverse-correlation invariants before running.
