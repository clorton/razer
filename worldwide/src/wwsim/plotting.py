"""Validation plots for the worldwide pipeline.

Liberal, illustrative figures: a global admin-2 population choropleth, the cross-border air
network drawn over the world, the airport->node assignment, and the air-edge weight
distribution. All write PNGs into ``output/plots/`` and return the Matplotlib figure.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: render to files, no display needed
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

from .logging import logger  # noqa: E402


def plot_population_choropleth(nodes_gdf, out_path: Path) -> plt.Figure:
    """World admin-2 population choropleth (log10 scale).

    Args:
        nodes_gdf: Global node GeoDataFrame with ``population`` and polygon geometry.
        out_path: PNG output path.

    Returns:
        The Matplotlib figure.
    """
    fig, ax = plt.subplots(figsize=(20, 10))
    g = nodes_gdf.copy()
    g["log_pop"] = np.log10(g["population"].clip(lower=1))
    g.plot(column="log_pop", ax=ax, cmap="inferno", linewidth=0, legend=True,
           legend_kwds={"label": "log10(population)", "shrink": 0.5})
    ax.set_title(f"Global admin-2 population ({g['iso3'].nunique()} countries, "
                 f"{len(g):,} units, {g['population'].sum()/1e9:.2f}B people)")
    ax.set_axis_off()
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    logger.info("plot: wrote %s", out_path)
    return fig


def plot_air_network(nodes_gdf, air_edges: pd.DataFrame, out_path: Path,
                     max_edges: int = 4000) -> plt.Figure:
    """Cross-border air network drawn over a faint world admin-2 background.

    Args:
        nodes_gdf: Global node GeoDataFrame (provides ``global_nodeid, lon, lat``).
        air_edges: Air network edges with ``src/dst_global_nodeid, weight``.
        out_path: PNG output path.
        max_edges: Draw at most this many heaviest edges (keeps the figure legible).

    Returns:
        The Matplotlib figure.
    """
    coords = nodes_gdf.set_index("global_nodeid")[["lon", "lat"]]
    edges = air_edges.sort_values("weight", ascending=False).head(max_edges)

    fig, ax = plt.subplots(figsize=(20, 10))
    nodes_gdf.plot(ax=ax, color="0.92", edgecolor="0.8", linewidth=0.1)

    segs, widths = [], []
    wmax = edges["weight"].max() or 1.0
    for _, e in edges.iterrows():
        try:
            s = coords.loc[e["src_global_nodeid"]]
            d = coords.loc[e["dst_global_nodeid"]]
        except KeyError:
            continue
        segs.append([(s["lon"], s["lat"]), (d["lon"], d["lat"])])
        widths.append(0.1 + 2.5 * (e["weight"] / wmax))
    lc = LineCollection(segs, colors="crimson", linewidths=widths, alpha=0.35)
    ax.add_collection(lc)
    ax.set_title(f"Global cross-border air network (top {len(segs):,} of {len(air_edges):,} edges)")
    ax.set_axis_off()
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    logger.info("plot: wrote %s", out_path)
    return fig


def plot_airport_assignment(nodes_gdf, assignment: pd.DataFrame, airports: pd.DataFrame,
                            out_path: Path) -> plt.Figure:
    """Airports plotted over the world, colored by assignment method.

    Args:
        nodes_gdf: Global node GeoDataFrame (world background).
        assignment: Output of :func:`wwsim.airports_assign.assign_airports_to_nodes`.
        airports: Airport master (provides ``iata, lon, lat``).
        out_path: PNG output path.

    Returns:
        The Matplotlib figure.
    """
    ap = airports.reset_index() if "iata" not in airports.columns else airports
    merged = assignment.merge(ap[["iata", "lon", "lat"]], on="iata", how="left")

    fig, ax = plt.subplots(figsize=(20, 10))
    nodes_gdf.plot(ax=ax, color="0.95", edgecolor="0.85", linewidth=0.1)
    for method, color in (("within", "navy"), ("nearest", "orange")):
        sub = merged[merged["assign_method"] == method]
        ax.scatter(sub["lon"], sub["lat"], s=2, c=color, label=f"{method} ({len(sub)})", alpha=0.6)
    ax.legend(markerscale=4, loc="lower left")
    ax.set_title(f"Airport -> admin-2 assignment ({len(merged):,} airports)")
    ax.set_axis_off()
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    logger.info("plot: wrote %s", out_path)
    return fig


def plot_edge_weight_distribution(air_edges: pd.DataFrame, out_path: Path) -> plt.Figure:
    """Histogram (log-log) of cross-border air-edge weights.

    Args:
        air_edges: Air network edges with ``weight``.
        out_path: PNG output path.

    Returns:
        The Matplotlib figure.
    """
    fig, ax = plt.subplots(figsize=(9, 6))
    w = air_edges["weight"].to_numpy()
    w = w[w > 0]
    ax.hist(w, bins=np.logspace(np.log10(w.min()), np.log10(w.max()), 50), color="steelblue")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("edge weight (passenger-proxy: seats offered)")
    ax.set_ylabel("number of admin-2 edges")
    ax.set_title("Cross-border air-edge weight distribution")
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    logger.info("plot: wrote %s", out_path)
    return fig


def plot_country_gravity(nodes_gdf, iso3: str, ids: np.ndarray, matrix: np.ndarray,
                         out_path: Path) -> plt.Figure:
    """One country's admin-2 gravity network drawn over its boundaries.

    Args:
        nodes_gdf: Global node GeoDataFrame.
        iso3: Country to draw.
        ids: Global node ids for the matrix rows/columns.
        matrix: Dense gravity matrix for the country.
        out_path: PNG output path.

    Returns:
        The Matplotlib figure.
    """
    country = nodes_gdf[nodes_gdf["iso3"] == iso3]
    coords = nodes_gdf.set_index("global_nodeid")[["lon", "lat"]]

    fig, ax = plt.subplots(figsize=(10, 10))
    country.plot(ax=ax, color="0.95", edgecolor="0.7", linewidth=0.3)

    # Draw the strongest few edges to keep it readable.
    flat = matrix.copy()
    thresh = np.quantile(flat[flat > 0], 0.95) if (flat > 0).any() else 0
    segs, widths = [], []
    wmax = flat.max() or 1.0
    rows, cols = np.where(flat >= thresh)
    for r, c in zip(rows, cols):
        s, d = coords.loc[ids[r]], coords.loc[ids[c]]
        segs.append([(s["lon"], s["lat"]), (d["lon"], d["lat"])])
        widths.append(0.1 + 2.0 * flat[r, c] / wmax)
    ax.add_collection(LineCollection(segs, colors="darkgreen", linewidths=widths, alpha=0.3))
    ax.scatter(country["lon"], country["lat"], s=6, c="black", zorder=3)
    ax.set_title(f"{iso3}: intra-country gravity network ({len(country)} admin-2 units)")
    ax.set_axis_off()
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    logger.info("plot: wrote %s", out_path)
    return fig
