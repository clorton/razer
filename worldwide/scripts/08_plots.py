"""Generate validation plots into output/plots/.

Example:
    python scripts/08_plots.py
    python scripts/08_plots.py --country NGA
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim import plotting  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.gravity import build_country_matrix  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--country", default="NGA", help="ISO3 for the per-country gravity plot")
    args = p.parse_args(argv)

    cfg = load_config()
    cfg.ensure_dirs()
    nodes = load_global_nodes(cfg)

    plotting.plot_population_choropleth(nodes, cfg.plots_dir / "global_population_choropleth.png")

    air_path = cfg.networks_dir / "global_air_network.parquet"
    if air_path.exists():
        air = pd.read_parquet(air_path)
        plotting.plot_air_network(nodes, air, cfg.plots_dir / "global_air_network.png")
        plotting.plot_edge_weight_distribution(air, cfg.plots_dir / "air_edge_weights.png")

    assign_path = cfg.interim_dir / "airport_node_assignment.parquet"
    airports_path = cfg.interim_dir / "airports_master.parquet"
    if assign_path.exists() and airports_path.exists():
        plotting.plot_airport_assignment(
            nodes, pd.read_parquet(assign_path), pd.read_parquet(airports_path),
            cfg.plots_dir / "airport_assignment.png",
        )

    if (nodes["iso3"] == args.country).any():
        grp = nodes[nodes["iso3"] == args.country]
        ids, matrix = build_country_matrix(grp, cfg.gravity)
        plotting.plot_country_gravity(
            nodes, args.country, ids, matrix, cfg.plots_dir / f"{args.country}_gravity.png"
        )

    logger.info("plots: done -> %s", cfg.plots_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
