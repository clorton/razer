"""Animate the simulation as 'new cases' on the SAME scale/aggregation as the COVID movie.

Renders the sim's **daily new infections per 100k (7-day avg)**, aggregated by country, with
the *identical* color scale (vmax), colormap, norm, and country aggregation as
``scripts/13_animate_covid.py`` -- so the two mp4s are directly comparable side by side.

The sim records compartments, not the incidence flow, but in a plain SEIR (no waning) the
daily new infections per node are exactly the daily drop in S: ``incidence[t] = S[t]-S[t+1]``.

The shared color ceiling is the COVID series' 99th-percentile per-100k value (computed here
from the OWID data, identical to what the COVID movie used), so the sim is shown on the real
world's scale.

Example:
    python scripts/14_animate_sim_cases.py --history output/seir/history_world540_air.npz --seconds 30
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim import covid_reference as cov  # noqa: E402
from wwsim.choropleth_anim import country_vmax, render_country_movie  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history", type=Path, default=Path("output/seir/history_world540_air.npz"))
    p.add_argument("--seconds", type=float, default=30.0)
    p.add_argument("--fps", type=int, default=None)
    p.add_argument("--width", type=int, default=3000)
    p.add_argument("--cmap", default="inferno")
    p.add_argument("--covid-days", type=int, default=730, help="window used to set the shared vmax")
    p.add_argument("--vmax", type=float, default=None, help="override the shared color ceiling")
    p.add_argument("--out", type=Path, default=None)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    cfg = load_config()
    out = args.out or (cfg.output_dir / "seir" / f"anim_cases_{args.history.stem.replace('history_', '')}.mp4")

    nodes = load_global_nodes(cfg).sort_values("global_nodeid").reset_index(drop=True)
    iso3_order = sorted(nodes["iso3"].unique())
    order_idx = {iso: k for k, iso in enumerate(iso3_order)}
    node_country_idx = nodes["iso3"].map(order_idx).to_numpy()
    n_countries = len(iso3_order)
    # One-hot node -> country, to aggregate per-node arrays to per-country (sum over nodes).
    onehot = sp.csr_matrix((np.ones(len(nodes)), (np.arange(len(nodes)), node_country_idx)),
                           shape=(len(nodes), n_countries))

    # --- shared color ceiling from the real COVID data (same as the COVID movie) ---
    full_pop = nodes.groupby("iso3")["population"].sum().to_dict()
    covid_df = cov.load_owid(cfg)
    _, covid_vals = cov.country_value_matrix(covid_df, iso3_order, full_pop,
                                             n_days=args.covid_days, field="new_cases")
    vmax = args.vmax if args.vmax is not None else country_vmax(covid_vals)
    logger.info("sim-cases: shared color ceiling vmax=%.2f per-100k/day (from COVID p99)", vmax)

    # --- sim daily new infections per node = drop in S ---
    hist = np.load(args.history)
    S = hist["S"].astype(np.int64)                  # (nticks+1, n_nodes)
    incidence_node = np.clip(-np.diff(S, axis=0), 0, None).astype(np.float64)  # (nticks, n_nodes)

    # Aggregate to country, convert to per-100k, smooth 7-day -- matching the COVID metric.
    inc_country = (onehot.T @ incidence_node.T).T   # (nticks, n_countries)
    country_pop_sub = np.asarray(onehot.T @ hist["node_pop"]).ravel()  # subsampled country pop
    rate = inc_country / np.maximum(country_pop_sub, 1)[None, :] * 100_000
    rate = pd.DataFrame(rate).rolling(7, min_periods=1).mean().to_numpy().astype(np.float32)

    frame_labels = [f"simulation — daily new infections — day {t + 1}" for t in range(rate.shape[0])]
    fps = args.fps or max(1, round(rate.shape[0] / args.seconds))

    render_country_movie(
        nodes, node_country_idx, rate, frame_labels, out,
        vmax=vmax, colorbar_label="daily new infections per 100k (7-day avg)",
        cmap=args.cmap, gamma=0.4, width=args.width, fps=fps,
    )
    logger.info("sim-cases: sim peak country rate=%.0f per-100k/day (ceiling %.0f) -> saturates above ceiling",
                float(np.nanmax(rate)), vmax)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
