"""Plots for a worldwide SEIR run: global epi curve, attack-rate map, country curves."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..logging import logger  # noqa: E402


def plot_epi_curve(totals: pd.DataFrame, out_path: Path, title: str = "Global SEIR") -> None:
    """Global S/E/I/R stacks and daily incidence over time.

    Args:
        totals: Output of :meth:`wwsim.abm.model.WorldSEIR.totals`.
        out_path: PNG path.
        title: Figure title.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    for st, c in (("S", "tab:blue"), ("E", "tab:purple"), ("I", "tab:orange"), ("R", "tab:green")):
        ax1.plot(totals.index, totals[st], label=st, color=c, lw=2)
    ax1.set_ylabel("people")
    ax1.legend(loc="center right")
    ax1.set_title(title)
    ax2.plot(totals.index, totals["incidence"], color="crimson", lw=1.5)
    ax2.set_ylabel("daily new infections")
    ax2.set_xlabel("day")
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    plt.close(fig)
    logger.info("abm.plots: wrote %s", out_path)


def plot_attack_choropleth(nodes_gdf, attack_rate: np.ndarray, out_path: Path) -> None:
    """World choropleth of the per-node final attack rate (R_final / population).

    Args:
        nodes_gdf: Global node GeoDataFrame (geometry; rows aligned to ``global_nodeid``).
        attack_rate: Per-node attack rate in [0, 1], length == number of nodes.
        out_path: PNG path.
    """
    g = nodes_gdf.copy()
    g["attack_rate"] = attack_rate
    fig, ax = plt.subplots(figsize=(20, 10))
    g.plot(column="attack_rate", ax=ax, cmap="inferno", vmin=0, vmax=1, linewidth=0,
           legend=True, legend_kwds={"label": "final attack rate", "shrink": 0.5})
    ax.set_title(f"Final attack rate by admin-2 (mean {np.average(attack_rate, weights=g['population']):.2f})")
    ax.set_axis_off()
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    plt.close(fig)
    logger.info("abm.plots: wrote %s", out_path)


def plot_country_curves(model, isos: list[str], out_path: Path) -> None:
    """Infectious-prevalence curves for selected countries (shows spatial/air spread timing).

    Args:
        model: A run :class:`wwsim.abm.model.WorldSEIR`.
        isos: ISO3 codes to plot.
        out_path: PNG path.
    """
    iso_col = model.scenario["iso3"].to_numpy()
    fig, ax = plt.subplots(figsize=(12, 7))
    for iso in isos:
        cols = np.where(iso_col == iso)[0]
        if len(cols) == 0:
            continue
        ax.plot(model.nodes.I[:, cols].sum(axis=1), label=iso, lw=2)
    ax.set_xlabel("day")
    ax.set_ylabel("infectious")
    ax.set_title("Infectious prevalence by country (arrival timing)")
    ax.legend(ncol=2, fontsize=8)
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    plt.close(fig)
    logger.info("abm.plots: wrote %s", out_path)
