"""Validated CSV loading reused from the amOT MDS package."""
import csv
import re
import numpy as np
import pandas as pd

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
    return pd.DataFrame(values, columns=keys), genes, len(df)

