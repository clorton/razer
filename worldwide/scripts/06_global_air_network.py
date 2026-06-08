"""Build the global cross-border air-travel network over admin-2 nodes.

Reads airport edges (step 03) and the airport->node assignment (step 04), and writes
``output/networks/global_air_network.{parquet,npz}``.

The ``--top-n`` option rebuilds the network for one or more airport cut-offs (busiest N by
passenger-volume proxy). A single value writes the canonical ``global_air_network.*``;
several values each write a suffixed ``global_air_network_top<N>.*`` so they coexist. This
step is cheap and independent of the node table and per-country gravity matrices.

Examples:
    python scripts/06_global_air_network.py                 # default top-N from config (1000)
    python scripts/06_global_air_network.py --top-n 500      # busiest 500 airports
    python scripts/06_global_air_network.py --top-n 1000 500 250 100   # all four variants
    python scripts/06_global_air_network.py --top-n 0        # 0 / "all" = no airport cap
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.air_network import build_air_network, save_air_network  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--top-n", nargs="+", type=int, default=None,
                   help="airport cut-off(s); 0 means no cap (keep all airports)")
    args = p.parse_args(argv)

    cfg = load_config()
    nodes = load_global_nodes(cfg)
    airport_edges = pd.read_parquet(cfg.interim_dir / "airport_edges.parquet")
    assignment = pd.read_parquet(cfg.interim_dir / "airport_node_assignment.parquet")

    # Default: one network at the config's top-N, canonical filename.
    top_values = args.top_n if args.top_n is not None else [cfg.air.top_n_airports]
    multiple = len(top_values) > 1

    for n in top_values:
        cfg.air.top_n_airports = None if (n is None or n == 0) else n
        edges = build_air_network(airport_edges, assignment, cfg)
        label = "all" if cfg.air.top_n_airports is None else cfg.air.top_n_airports
        stem = f"global_air_network_top{label}" if multiple else "global_air_network"
        save_air_network(edges, cfg, n_nodes=len(nodes), stem=stem)
        logger.info(
            "air_network[top=%s]: %d cross-border admin-2 edges, %d gateway nodes, weight=%s",
            label, len(edges), edges["src_global_nodeid"].nunique(), f"{edges['weight'].sum():,.0f}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
