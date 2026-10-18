from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt


def finish_figure(fig, output=None, dpi: int = 300, close: bool = False):
    fig.tight_layout()
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(fig)
    return fig


def save_formats(fig, base_path, formats=("png", "pdf"), dpi: int = 300):
    """Save an existing plot in raster and vector formats without recalculation."""
    base = Path(base_path).with_suffix("")
    return [finish_figure(fig, base.with_suffix("."+fmt), dpi) for fmt in formats]
