"""Plot five pooled right-amOT mean profiles from prepared cell-by-gene CSVs."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import sys
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata
from scipy.spatial.distance import pdist, squareform
import sklearn
from sklearn.manifold import MDS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from plot_helpers.embedding import embedding_plot

GROUPS = [('tD1', 'c61'), ('D2', 'c62'), ('Sema5a-D1', 'c63'),
          ('pDwarf Cell', 'c264'), ('ICj', 'c262')]
METADATA = {'cell_label', 'cluster_alias', 'brain_section_label', 'x', 'y', 'z',
            'x_ccf', 'y_ccf', 'z_ccf', 'parcellation_structure',
            'parcellation_substructure', 'section', 'region_id'}


def gene_key(name):
    """Reconcile MATLAB-valid column names, requiring a bijection downstream."""
    key = re.sub(r'[^A-Za-z0-9_]', '_', name)
    return 'x' + key if key[0].isdigit() else key


def load_profile(path):
    """Validate the pre-partitioned right-amOT input and average supplied values."""
    with path.open(encoding='utf-8-sig', newline='') as stream:
        header = next(csv.reader(stream))
    if len(header) != len(set(header)):
        raise ValueError(f'Duplicate column names: {path.name}')
    df = pd.read_csv(path, dtype={'cell_label': str})
    for column, value in [('partition_domain', 'amOT'), ('partition_hemisphere', 'R')]:
        if column not in df or not df[column].eq(value).all():
            raise ValueError(f'{path.name}: every row must have {column}={value}')
    if df.empty or 'z_ccf' not in df or not (pd.to_numeric(df.z_ccf) < 5.7).all():
        raise ValueError(f'{path.name}: expected nonempty right-amOT cells, z_ccf < 5.7')
    if 'partition_source_row' not in df or df.partition_source_row.isna().any() or not df.partition_source_row.is_unique:
        raise ValueError(f'{path.name}: partition_source_row must be present and unique')
    genes = [g for g in header if g not in METADATA and not g.startswith('partition_')]
    keys = [gene_key(g) for g in genes]
    if len(genes) != 8460 or len(set(keys)) != 8460:
        raise ValueError(f'{path.name}: expected 8460 unambiguous gene columns, found {len(genes)}')
    values = df[genes].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f'{path.name}: expression contains missing/infinite values')
    return pd.Series(values.mean(axis=0), index=keys), genes, len(df)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True, help='Directory containing the five *_OT_R_amOT_v3.csv files')
    parser.add_argument('--output-dir', type=Path, default=Path('results'))
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    profiles, headers, sources = [], [], []
    for name, prefix in GROUPS:
        path = args.input_dir / f'{prefix}_OT_R_amOT_v3.csv'
        profile, genes, n = load_profile(path)
        profiles.append(profile)
        headers.append(genes)
        sources.append({'group': name, 'file': path.name, 'n_cells': n,
                        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        print(f'{name}: {n} cells, {len(genes)} genes', flush=True)
    genes = headers[1]  # D2 export retains the canonical gene symbols.
    keys = [gene_key(g) for g in genes]
    if any(set(p.index) != set(keys) for p in profiles):
        raise ValueError('The five gene sets do not match exactly after MATLAB-name reconciliation')
    x = np.stack([p.loc[keys].to_numpy() for p in profiles])
    names = [name for name, _ in GROUPS]
    pd.DataFrame(x.T, index=pd.Index(genes, name='gene'), columns=names).to_csv(out / 'mean_profiles.csv')
    mapping = []
    canonical = dict(zip(keys, genes))
    for name, columns in zip(names, headers):
        mapping.extend({'group': name, 'csv_column': g, 'gene_symbol': canonical[gene_key(g)]} for g in columns)
    pd.DataFrame(mapping).to_csv(out / 'gene_name_mapping.csv', index=False)
    ranks = rankdata(x, axis=1)
    ranks -= ranks.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranks, axis=1, keepdims=True)
    if (norms == 0).any():
        raise ValueError('A profile is constant across all genes')
    unit = ranks / norms
    rho = unit @ unit.T
    distance = np.sqrt(np.maximum(0, 2 * (1 - rho)))
    np.fill_diagonal(distance, 0)
    xy = MDS(n_components=2, metric=True, dissimilarity='precomputed', random_state=0,
             n_init=32, max_iter=3000, eps=1e-9, normalized_stress=False).fit_transform(distance)
    projected = squareform(pdist(xy))
    upper = np.triu_indices(5, 1)
    stress = float(np.sqrt(np.sum((projected[upper] - distance[upper])**2) / np.sum(distance[upper]**2)))
    for filename, matrix in [('spearman_correlations', rho), ('correlation_distances', distance)]:
        pd.DataFrame(matrix, index=names, columns=names).to_csv(out / f'{filename}.csv')
    pd.DataFrame(xy, index=names, columns=['MDS1', 'MDS2']).to_csv(out / 'coordinates.csv')
    plt.rcParams.update({'font.family': ['Arial', 'DejaVu Sans'], 'font.size': 14, 'svg.fonttype': 'none'})
    fig = embedding_plot(xy)
    fig.set_size_inches(5, 5)
    ax = fig.axes[0]
    ax.collections[0].remove()
    colors = ['#d45143', '#347ba5', '#338a77', '#b45a82', '#8267a5']
    offsets = [(0, 16), (-10, 10), (0, -17), (0, -18), (0, 17)]
    for point, name, color, offset in zip(xy, names, colors, offsets):
        ax.scatter(*point, s=110, c=color, edgecolors='white', linewidths=1, zorder=3)
        ax.annotate(name, point, xytext=offset, textcoords='offset points', ha='center', va='center', fontsize=16, color=color)
    # Original figure limits, expanded automatically for other validated inputs.
    ax.set_xlim(min(-.38, xy[:, 0].min() - .16), max(.38, xy[:, 0].max() + .16))
    ax.set_ylim(min(-.31, xy[:, 1].min() - .09), max(.45, xy[:, 1].max() + .07))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('MDS 1', fontsize=16)
    ax.set_ylabel('MDS 2', fontsize=16)
    available = {f.name for f in font_manager.fontManager.ttflist}
    title_font = next((f for f in ['Meiryo', 'Noto Sans CJK JP', 'Yu Gothic'] if f in available), 'DejaVu Sans')
    # Use a visually equivalent middle dot if no Japanese font is installed.
    title = 'amOT・R' if title_font != 'DejaVu Sans' else 'amOT·R'
    ax.set_title(title, fontsize=18, fontfamily=title_font, pad=12)
    ax.set_xticks([-.3, 0, .3])
    ax.set_yticks([-.2, 0, .2, .4])
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.14)
    fig.subplots_adjust(left=.19, right=.98, bottom=.15, top=.90)
    for ext in ['png', 'svg', 'pdf']:
        fig.savefig(out / f'amOT_R_MDS.{ext}', dpi=400)
    plt.close(fig)
    manifest = {'sources': sources, 'genes': len(genes), 'summary': 'mean of supplied log2 values; cells pooled across sections',
                'distance': 'sqrt(2*(1-Spearman rho))', 'method': 'metric MDS', 'seed': 0,
                'n_init': 32, 'max_iter': 3000, 'eps': 1e-9, 'normalized_distance_stress': stress,
                'notes': 'One point per group, not per cell or donor. Descriptive plot; no significance test. 2D distances approximate original distances. No imputation, normalization, or rounding performed.',
                'versions': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__, 'scikit-learn': sklearn.__version__, 'matplotlib': matplotlib.__version__}}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(f'Saved to {out.resolve()} | normalized distance stress = {stress:.6f}')


if __name__ == '__main__':
    main()
