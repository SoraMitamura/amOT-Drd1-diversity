# amOT Drd1 diversity

Analysis accompanying **Four populations of Drd1-expressing neurons with distinct transcriptomic profiles in the anteromedial olfactory tubercle**.

**Sora Mitamura**  
Division of Brain Structure and Function, Faculty of Medical Sciences, University of Fukui, Fukui 910-1193, Japan

## Scope and results

Transcriptomic similarity to typical D1 MSNs (tD1) in **right amOT v3**, pooling available sections. Comparators: D2, Sema5a-D1, pDwarf Cell (c264), ICj (c262). This repository contains downstream code, archived results, calculation records and figures. Input cell-by-gene matrices and the manuscript are not distributed here.

- [All-gene ranking](results/all_genes/similarity_ranking.csv)
- [Six tD1-reference contrasts](results/all_genes/test_results.csv)
- [Eight marker-excluded GO categories and figures](results/GO)
- [All 48 GO tests](results/GO/all_48_tests_marker_excluded.csv)
- [Before/after marker exclusion](results/GO/before_after_correlations.csv)
- [All-gene MDS code and outputs](mds)

![All-gene contrasts](results/all_genes/all_six_correlation_tests.png)

Across 8,460 analyzed genes, observed similarity to tD1 is D2 (0.9821), Sema5a-D1 (0.9761), pDwarf Cell (0.9686), then ICj (0.8643). All six raw p-values reach the Monte Carlo resolution floor of 0.0001, not a true probability of zero.

## Reproduce

```sh
python -m pip install -r requirements.txt
python all_genes/analysis.py --input-dir data --output-dir rerun/all_genes
python run_go.py --input-dir data --output-dir rerun/GO
```

Supply five original right-amOT CSVs in `data/`: `c61_OT_R_amOT_v3.csv` (912 cells), `c62_OT_R_amOT_v3.csv` (501), `c63_OT_R_amOT_v3.csv` (101), `c264_OT_R_amOT_v3.csv` (352), and `c262_OT_R_amOT_v3.csv` (468). These are not section-specific coronal exports. Input SHA256 hashes are in each result directory's `input_provenance.csv`. Loaders enforce right-amOT metadata, unique source rows, finite expression and matching gene panels with the existing MATLAB gene-name normalization. See [full input/method details](all_genes/README.md) and [MDS CLI instructions](mds/README.md).

## Methods and interpretation

Population profiles are means of supplied log2 expression across cells. Spearman correlations are computed across genes. The six contrasts are rho(tD1,A) minus rho(tD1,B). Whole cells are independently bootstrapped within populations, reusing the same sampled tD1 profile for all comparisons: 9,999 replicates, seed 20261017, batch size 64. Outputs include percentile 95% intervals and two-sided centered-bootstrap p-values with a plus-one correction, raw and six-test Holm p-values, and BH correction across the 48 GO tests. Saved tail indicators preserve the exact calculation despite CSV serialization.

MDS uses sqrt(2*(1-rho)), metric MDS, seed 0 and 32 initializations. MDS is descriptive; tests compare correlations, not two-dimensional positions.

GO definitions come from `m5.go.bp.v2026.1.Mm.symbols.gmt`, SHA256 `fd720ba92a3131a1596409f4b642a1a5950392fa51b45055b810215082e6f0c1`. Frozen selections retain category URLs and membership audits. This is category-restricted similarity, not enrichment analysis. The marker exclusion panel is **Drd1, Tac1, Pdyn, Drd2, Adora2a, Penk, Gpr6, Gpr52, Sp9**.

| Category | Genes after exclusion | Removed |
|---|---:|---|
| Neuron projection guidance | 196 | None |
| Synapse organization | 510 | Drd1, Drd2 |
| Dendrite development | 204 | None |
| Synaptic vesicle cycle | 156 | None |
| Neuronal action potential | 41 | Drd1 |
| Regulation of postsynaptic membrane potential | 144 | Adora2a |
| Neuropeptide signaling pathway | 86 | Pdyn, Penk, Tac1 |
| Gamma-aminobutyric acid signaling pathway | 33 | None |

These are exploratory comparisons conditional on **one spatial specimen's fixed imputed matrix**. Cells are not independent animal replicates. Cell exchangeability is assumed; spatial dependence and imputation uncertainty are not propagated. Small p-values are not animal-level evidence. No section stratification, new expression threshold, normalization or rounding is applied. pDwarf and ICj are operational annotation assignments, not newly established morphological identities. The 8,460 genes are the analyzed panel, not the entire mouse genome.

All 54 saved contrasts were checked against bootstrap intervals and tail counts; see [validation](validation.txt). Absolute workstation paths in JSON were replaced with source basenames, preserving input and GMT hashes.
