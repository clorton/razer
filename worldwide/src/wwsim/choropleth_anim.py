"""Shared country-aggregated choropleth movie renderer.

Both the real-COVID animation and the sim-as-cases animation use this so they share the
*exact* same rendering: same rasterization, color norm, colormap, and (crucially) the same
``vmax`` -- which is what makes the two movies directly comparable.

Approach (cheap per frame): rasterize the admin-2 polygons to a node-index grid once, then
each frame color every node by **its country's** value -- a country choropleth with no
polygon dissolve.
"""

from __future__ import annotations

from pathlib import Path

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import rasterio.features  # noqa: E402
from matplotlib.colors import PowerNorm  # noqa: E402
from rasterio.transform import from_bounds  # noqa: E402

from .logging import logger  # noqa: E402


def rasterize_nodes(nodes_gdf, width: int):
    """Burn admin-2 polygons into a grid of node ids (-1 = background), once.

    Args:
        nodes_gdf: Global node GeoDataFrame (row k == global_nodeid k).
        width: Output width in pixels; height follows the data aspect ratio.

    Returns:
        Tuple ``(idx int32 grid, bounds (minx, miny, maxx, maxy))``.
    """
    minx, miny, maxx, maxy = nodes_gdf.total_bounds
    height = max(1, int(round(width * (maxy - miny) / (maxx - minx))))
    transform = from_bounds(minx, miny, maxx, maxy, width, height)
    shapes = ((g, int(i)) for g, i in zip(nodes_gdf.geometry, nodes_gdf["global_nodeid"]))
    logger.info("choropleth_anim: rasterizing %d polygons to %dx%d...", len(nodes_gdf), width, height)
    idx = rasterio.features.rasterize(shapes, out_shape=(height, width), transform=transform,
                                      fill=-1, all_touched=True, dtype=np.int32)
    return idx, (minx, miny, maxx, maxy)


def country_vmax(values: np.ndarray, quantile: float = 0.99) -> float:
    """The color-scale ceiling: a high quantile of positive values (outlier-robust).

    Args:
        values: ``(T, n_countries)`` value matrix.
        quantile: Quantile of positive values to use.

    Returns:
        The vmax, or 1.0 if there are no positive values.
    """
    pos = values[np.isfinite(values) & (values > 0)]
    return float(np.quantile(pos, quantile)) if pos.size else 1.0


def render_country_movie(
    nodes_gdf,
    node_country_idx: np.ndarray,
    values: np.ndarray,
    frame_labels: list[str],
    out: Path,
    *,
    vmax: float,
    colorbar_label: str,
    cmap: str = "inferno",
    gamma: float = 0.4,
    width: int = 3000,
    fps: int = 18,
) -> float:
    """Render a country-aggregated choropleth movie to ``out`` (mp4).

    Args:
        nodes_gdf: Global node GeoDataFrame (geometry + global_nodeid, sorted by id).
        node_country_idx: Per-node column index into ``values`` (node -> its country).
        values: ``(T, n_countries)`` per-frame country values (NaN = no data).
        frame_labels: Per-frame title text (length T).
        out: Output mp4 path.
        vmax: Color-scale ceiling (shared across movies for comparability).
        colorbar_label: Colorbar label.
        cmap: Matplotlib colormap name.
        gamma: PowerNorm gamma (<1 lifts the low-value frontier).
        width: Rasterization width in pixels.
        fps: Frames per second.

    Returns:
        The ``vmax`` used (echoed for logging).
    """
    idx, (minx, miny, maxx, maxy) = rasterize_nodes(nodes_gdf, width)
    bg = idx < 0
    safe = np.where(bg, 0, idx)

    norm = PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("0.85")  # background and no-data shown light gray

    def grid(t: int):
        node_vals = values[t][node_country_idx]
        g = node_vals[safe]
        return np.ma.masked_array(g, mask=bg | ~np.isfinite(g))

    fig, ax = plt.subplots(figsize=(16, 9), dpi=120)
    im = ax.imshow(grid(0), cmap=cm, norm=norm, extent=(minx, maxx, miny, maxy),
                   origin="upper", interpolation="nearest")
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01).set_label(colorbar_label)
    title = ax.set_title(frame_labels[0], fontsize=14)
    fig.tight_layout()

    logger.info("choropleth_anim: writing %d frames at %d fps (vmax=%.3g) -> %s",
                len(values), fps, vmax, out)
    with imageio.get_writer(out, fps=fps, codec="libx264", quality=8, macro_block_size=16) as writer:
        for t in range(len(values)):
            im.set_data(grid(t))
            title.set_text(frame_labels[t])
            fig.canvas.draw()
            writer.append_data(np.asarray(fig.canvas.buffer_rgba())[..., :3])
    plt.close(fig)
    logger.info("choropleth_anim: done (%.1fs at %d fps) -> %s", len(values) / fps, fps, out)
    return vmax
