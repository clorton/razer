"""Animate an SEIR run as a choropleth movie (mp4) from a saved per-node history.

Reads ``output/seir/history_<tag>.npz`` (written by ``10_run_seir.py --save-history``) and
the global node geometry, and renders one frame per tick into an mp4.

Speed trick: drawing 45k polygons per frame would take hours, so we **rasterize the polygons
to a node-index grid once** and then, each frame, recolor that grid by indexing the chosen
field -- effectively free per frame.

The default field is **infectious prevalence** (I / population), which shows the wave moving;
``attack`` (R / population) shows the cumulative footprint growing.

Example:
    python scripts/12_animate_choropleth.py --history output/seir/history_world540_air.npz \
        --field I --seconds 30
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import rasterio.features  # noqa: E402
from matplotlib.colors import PowerNorm  # noqa: E402
from rasterio.transform import from_bounds  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history", type=Path, required=True, help="history_<tag>.npz from step 10")
    p.add_argument("--field", choices=["I", "E", "EI", "attack"], default="I",
                   help="I/E/EI = prevalence (rate); attack = cumulative R/pop")
    p.add_argument("--seconds", type=float, default=30.0, help="target movie duration")
    p.add_argument("--fps", type=int, default=None, help="override frames/sec")
    p.add_argument("--width", type=int, default=3000, help="rasterization width in pixels")
    p.add_argument("--cmap", default="magma")
    p.add_argument("--out", type=Path, default=None)
    return p.parse_args(argv)


def field_series(hist, field: str) -> tuple[np.ndarray, str]:
    """Return the per-tick per-node field array and a label.

    Args:
        hist: Loaded npz with S/E/I/R and node_pop.
        field: One of I, E, EI, attack.

    Returns:
        Tuple of ``(values[nticks+1, n_nodes] float32, label)`` -- a rate in [0, ~].
    """
    pop = np.maximum(hist["node_pop"], 1.0)
    if field == "attack":
        return (hist["R"] / pop).astype(np.float32), "attack rate (R / population)"
    counts = hist["I"] if field == "I" else hist["E"] if field == "E" else (hist["E"] + hist["I"])
    label = {"I": "infectious", "E": "exposed", "EI": "exposed+infectious"}[field]
    return (counts / pop).astype(np.float32), f"{label} prevalence (/ population)"


def rasterize_nodes(nodes_gdf, width: int) -> tuple[np.ndarray, tuple]:
    """Burn admin-2 polygons into a grid of node ids (-1 = background), once.

    Args:
        nodes_gdf: Global node GeoDataFrame (row k == global_nodeid k).
        width: Output width in pixels; height follows the data aspect ratio.

    Returns:
        Tuple of (``idx`` int32 grid, bounds tuple ``(minx, miny, maxx, maxy)``).
    """
    minx, miny, maxx, maxy = nodes_gdf.total_bounds
    height = max(1, int(round(width * (maxy - miny) / (maxx - minx))))
    transform = from_bounds(minx, miny, maxx, maxy, width, height)
    shapes = ((geom, int(nid)) for geom, nid in zip(nodes_gdf.geometry, nodes_gdf["global_nodeid"]))
    logger.info("animate: rasterizing %d polygons to %dx%d (one time)...", len(nodes_gdf), width, height)
    idx = rasterio.features.rasterize(
        shapes, out_shape=(height, width), transform=transform, fill=-1,
        all_touched=True, dtype=np.int32,
    )
    return idx, (minx, miny, maxx, maxy)


def main(argv=None) -> int:
    args = parse_args(argv)
    cfg = load_config()
    hist = np.load(args.history)
    nodes = load_global_nodes(cfg).sort_values("global_nodeid").reset_index(drop=True)

    values, label = field_series(hist, args.field)
    n_frames = values.shape[0]
    fps = args.fps or max(1, round(n_frames / args.seconds))
    out = args.out or (cfg.output_dir / "seir" / f"anim_{args.history.stem.replace('history_', '')}_{args.field}.mp4")

    idx, bounds = rasterize_nodes(nodes, args.width)
    bg = idx < 0
    safe = np.where(bg, 0, idx)

    # Fixed color scale across the whole movie so frames are comparable.
    positive = values[values > 0]
    vmax = float(np.quantile(positive, 0.99)) if positive.size else 1.0
    norm = PowerNorm(gamma=0.4, vmin=0.0, vmax=vmax)  # gamma<1 lifts the low-prevalence frontier
    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad("0.85")  # background land/sea shown light gray

    minx, miny, maxx, maxy = bounds
    fig, ax = plt.subplots(figsize=(16, 9), dpi=120)
    grid0 = np.ma.masked_where(bg, values[0][safe])
    im = ax.imshow(grid0, cmap=cmap, norm=norm, extent=(minx, maxx, miny, maxy), origin="upper", interpolation="nearest")
    ax.set_axis_off()
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label(label)
    title = ax.set_title(f"{label} — day 0", fontsize=14)
    fig.tight_layout()

    logger.info("animate: writing %d frames at %d fps -> %s", n_frames, fps, out)
    with imageio.get_writer(out, fps=fps, codec="libx264", quality=8, macro_block_size=16) as writer:
        for t in range(n_frames):
            im.set_data(np.ma.masked_where(bg, values[t][safe]))
            title.set_text(f"{label} — day {t}")
            fig.canvas.draw()
            frame = np.asarray(fig.canvas.buffer_rgba())[..., :3]
            writer.append_data(frame)
    plt.close(fig)
    logger.info("animate: done (%.1fs at %d fps) -> %s", n_frames / fps, fps, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
