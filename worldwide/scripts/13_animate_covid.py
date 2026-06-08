"""Animate real COVID-19 data (Our World in Data) as a country choropleth movie.

The real series are country-level, so this is a per-country choropleth. To avoid dissolving
45k admin-2 polygons, we reuse the admin-2 node raster and color each node by **its
country's** value -- which renders as a country choropleth (internal admin-2 borders vanish
because neighbours within a country share a colour).

Field: daily new cases per 100k (7-day average), the standard COVID choropleth metric. The
default window is the first 730 days from the dataset start (2020-01-05 → 2022-01-03:
ancestral → Alpha → Delta → early Omicron).

Example:
    python scripts/13_animate_covid.py --days 730 --seconds 30
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
import pandas as pd  # noqa: E402
import rasterio.features  # noqa: E402
from matplotlib.colors import PowerNorm  # noqa: E402
from rasterio.transform import from_bounds  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim import covid_reference as cov  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--days", type=int, default=730)
    p.add_argument("--field", choices=["new_cases", "new_deaths"], default="new_cases")
    p.add_argument("--per", type=float, default=100_000, help="per-capita denominator")
    p.add_argument("--seconds", type=float, default=30.0)
    p.add_argument("--fps", type=int, default=None)
    p.add_argument("--width", type=int, default=3000)
    p.add_argument("--cmap", default="inferno")
    p.add_argument("--out", type=Path, default=None)
    return p.parse_args(argv)


def rasterize_nodes(nodes_gdf, width: int):
    """Burn admin-2 polygons into a node-id grid (-1 background); see scripts/12."""
    minx, miny, maxx, maxy = nodes_gdf.total_bounds
    height = max(1, int(round(width * (maxy - miny) / (maxx - minx))))
    transform = from_bounds(minx, miny, maxx, maxy, width, height)
    shapes = ((g, int(i)) for g, i in zip(nodes_gdf.geometry, nodes_gdf["global_nodeid"]))
    logger.info("animate-covid: rasterizing %d polygons to %dx%d...", len(nodes_gdf), width, height)
    idx = rasterio.features.rasterize(shapes, out_shape=(height, width), transform=transform,
                                      fill=-1, all_touched=True, dtype=np.int32)
    return idx, (minx, miny, maxx, maxy)


def main(argv=None) -> int:
    args = parse_args(argv)
    cfg = load_config()
    out_dir = cfg.output_dir / "covid"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out or (out_dir / f"anim_covid_{args.days}d_{args.field}.mp4")

    nodes = load_global_nodes(cfg).sort_values("global_nodeid").reset_index(drop=True)
    iso3_order = sorted(nodes["iso3"].unique())
    populations = nodes.groupby("iso3")["population"].sum().to_dict()
    order_idx = {iso: k for k, iso in enumerate(iso3_order)}
    node_country = nodes["iso3"].map(order_idx).to_numpy()  # node -> column in iso3_order

    cov.download_owid(cfg)
    df = cov.load_owid(cfg)
    dates, values = cov.country_value_matrix(
        df, iso3_order, populations, n_days=args.days, field=args.field, per_capita_per=args.per
    )

    idx, bounds = rasterize_nodes(nodes, args.width)
    bg = idx < 0
    safe = np.where(bg, 0, idx)

    positive = values[np.isfinite(values) & (values > 0)]
    vmax = float(np.quantile(positive, 0.99)) if positive.size else 1.0
    norm = PowerNorm(gamma=0.4, vmin=0.0, vmax=vmax)
    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad("0.85")

    unit = "cases" if args.field == "new_cases" else "deaths"
    label = f"daily new COVID-19 {unit} per {int(args.per/1000)}k (7-day avg)"

    def frame_grid(t: int):
        node_vals = values[t][node_country]          # per-node = its country's value
        grid = node_vals[safe]
        return np.ma.masked_array(grid, mask=bg | ~np.isfinite(grid))

    minx, miny, maxx, maxy = bounds
    fig, ax = plt.subplots(figsize=(16, 9), dpi=120)
    im = ax.imshow(frame_grid(0), cmap=cmap, norm=norm, extent=(minx, maxx, miny, maxy),
                   origin="upper", interpolation="nearest")
    ax.set_axis_off()
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label(label)
    title = ax.set_title("", fontsize=14)
    fig.tight_layout()

    n_frames = args.days
    fps = args.fps or max(1, round(n_frames / args.seconds))
    logger.info("animate-covid: writing %d frames at %d fps -> %s", n_frames, fps, out)
    with imageio.get_writer(out, fps=fps, codec="libx264", quality=8, macro_block_size=16) as writer:
        for t in range(n_frames):
            im.set_data(frame_grid(t))
            title.set_text(f"{label}\n{pd.Timestamp(dates[t]).date()}")
            fig.canvas.draw()
            writer.append_data(np.asarray(fig.canvas.buffer_rgba())[..., :3])
    plt.close(fig)
    logger.info("animate-covid: done (%.1fs at %d fps) -> %s", n_frames / fps, fps, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
