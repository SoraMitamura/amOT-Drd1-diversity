# Transcriptomic similarity of five right-amOT populations

Generate the **amOT・R** metric MDS plot directly from five prepared cell-by-gene CSVs. No personal paths or external repository are required.

## Setup and run

Use Python 3.11 or 3.12. Install dependencies:

```bash
python -m pip install -r requirements.txt
python plot_amot_mds.py --input-dir data --output-dir results
```

Place these files in `data/` (or pass another directory):

| File | Display label |
|---|---|
| c61_OT_R_amOT_v3.csv | tD1 |
| c62_OT_R_amOT_v3.csv | D2 |
| c63_OT_R_amOT_v3.csv | Sema5a-D1 |
| c264_OT_R_amOT_v3.csv | pDwarf Cell |
| c262_OT_R_amOT_v3.csv | ICj |

The c264 label follows the population definition supplied for this analysis; the script does not infer cell identities from marker expression. Input data are not included in this repository.

## Input contract

Rows are cells and 8,460 columns are supplied imputed log2 gene expression values. Inputs must already be partitioned using the v3 procedure, with `partition_domain=amOT`, `partition_hemisphere=R`, `z_ccf<5.7`, and unique, nonmissing `partition_source_row` within each file. This script validates those fields; it does not perform anatomical partitioning.

Known metadata columns (`cell_label`, `cluster_alias`, `brain_section_label`, `x`, `y`, `z`, `x_ccf`, `y_ccf`, `z_ccf`, `parcellation_structure`, `parcellation_substructure`, `section`, `region_id`) and all `partition_` columns are excluded. Unknown extra columns cause validation to fail rather than silently entering the analysis.

MATLAB-style gene names are reconciled using punctuation-to-underscore conversion and an `x` prefix for numeric-leading symbols. The D2 file supplies canonical symbols. A one-to-one correspondence across all five files is required; missing genes, ambiguous names, duplicate headers, and nonfinite expression values stop execution. Rows are not deduplicated using `cell_label`.

## Methods

1. Pool all supplied right-amOT cells within each population, without section stratification.
2. Compute the arithmetic mean of the supplied log2 expression values for each gene. This is **not** log2 of mean linear expression. No transformation, normalization, additional rounding, or gene filtering is performed.
3. Compute Spearman correlations between five mean profiles over all 8,460 genes, with average ranks for ties.
4. Convert correlations to distances: `d = sqrt(2 * (1 - rho))`.
5. Apply metric MDS (two dimensions, seed 0, 32 initializations, 3,000 maximum iterations, tolerance 1e-9).

Each point is one population mean, not a cell or biological replicate. The visualization is descriptive and includes no statistical significance test. The two-dimensional positions approximate the original distances; rotation and axis signs have no biological meaning. The output manifest reports normalized distance stress, defined as `sqrt(sum((d_2D-d_original)^2)/sum(d_original^2))` over unique pairs. It is approximately 0.118 for the original dataset.

## Outputs

- `amOT_R_MDS.png` (400 dpi), `.svg`, `.pdf`
- `mean_profiles.csv`, `gene_name_mapping.csv`
- `spearman_correlations.csv`, `correlation_distances.csv`, `coordinates.csv`
- `manifest.json`: input hashes, cell counts, parameters, software versions, and stress

Arial is used when installed, with DejaVu Sans as fallback. The title uses a Japanese-capable font for `・` when available, otherwise the visually similar `·`. Input CSVs and generated results are excluded by `.gitignore`.

## Code provenance

The calculation and figure settings reproduce the previous amOT analysis. `plot_helpers/embedding.py` and `_common.py` reuse the existing NeuroResearchOS plotting helpers, bundled here to avoid a local-repository dependency. The main entry point is `plot_amot_mds.py`.
