"""Existing embeddings first; optional dimension reduction for bounded matrices."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ._common import finish_figure


def get_embedding(adata, *, method: str = "pca", existing_key: str | None = None,
                  matrix=None, seed: int = 0):
    """Return n_cells x 2 coordinates, using an h5ad obsm key when available.

    Matrix-based computation is opt-in: callers must pass a bounded matrix
    explicitly, with normalization and gene selection recorded in provenance.
    """
    key = existing_key or {"pca": "X_pca", "umap": "X_umap", "tsne": "X_tsne"}[method.lower()]
    if key in adata.obsm:
        return np.asarray(adata.obsm[key])[:, :2]
    if matrix is None:
        raise KeyError(f"Embedding {key} not present; supply an explicit bounded matrix")
    from sklearn.decomposition import PCA
    x = np.asarray(matrix, dtype=float)
    if x.ndim != 2 or min(x.shape) < 3 or not np.isfinite(x).all():
        raise ValueError("Embedding input must be finite cells x genes with >=3 rows/columns")
    if method.lower() == "pca":
        return PCA(n_components=2, random_state=seed).fit_transform(x)
    if method.lower() == "tsne":
        from sklearn.manifold import TSNE
        return TSNE(n_components=2, random_state=seed, perplexity=min(30, (len(x)-1)/3)).fit_transform(x)
    if method.lower() == "umap":
        import umap
        return umap.UMAP(n_components=2, random_state=seed).fit_transform(x)
    raise ValueError(f"Unknown embedding method: {method}")


def embedding_plot(coordinates, *, labels=None, title=None, output=None, dpi: int = 300):
    """Plot precomputed two-dimensional coordinates with optional group labels."""
    xy = np.asarray(coordinates)
    if xy.ndim != 2 or xy.shape[1] < 2:
        raise ValueError("Expected n x 2 coordinates")
    fig, ax = plt.subplots(figsize=(6, 5))
    if labels is None:
        ax.scatter(xy[:, 0], xy[:, 1], s=5, alpha=.6)
    else:
        labels = pd.Series(labels).reset_index(drop=True)
        if len(labels) != len(xy): raise ValueError("Label length differs from coordinates")
        for name in labels.unique():
            mask = (labels == name).to_numpy()
            ax.scatter(xy[mask, 0], xy[mask, 1], s=5, alpha=.6, label=str(name))
        ax.legend(markerscale=2)
    ax.set_xlabel("Component 1"); ax.set_ylabel("Component 2")
    if title: ax.set_title(title)
    return finish_figure(fig, output, dpi)
