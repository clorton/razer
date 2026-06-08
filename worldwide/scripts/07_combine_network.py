"""Combine intra-country gravity + inter-country air into one global migration matrix.

Rebuilds the per-country gravity matrices (fast) and the cross-border air network, then
sums them through the pluggable multi-modal combiner and writes a sparse ``.npz``.

``--top-n`` controls the airport cut-off for the air layer. A single value (or none) writes
the canonical ``global_combined_network.npz``; several values each write a suffixed
``global_combined_network_top<N>.npz`` so cut-offs coexist (mirrors step 06). When
``--top-n`` is given the air layer is rebuilt in-process from the step-03/04 interim files,
so step 06 need not have been run for that N.

``--gravity-scale`` / ``--air-scale`` make the two layers commensurate (gravity is raw
connection strength; air is the seat proxy). ``--rail-csv`` adds a rail mode (see
:mod:`wwsim.rail`).

Examples:
    python scripts/07_combine_network.py                       # canonical: config top-N air
    python scripts/07_combine_network.py --top-n 250            # canonical, top-250 air
    python scripts/07_combine_network.py --top-n 1000 500 250 100   # one combined file per N
    python scripts/07_combine_network.py --air-scale 1e9       # boost air relative to gravity
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.air_network import build_air_network  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.gravity import build_all_country_matrices  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.network import build_combined_network, save_combined_network  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--top-n", nargs="+", type=int, default=None,
                   help="airport cut-off(s) for the air layer; 0 means no cap (keep all)")
    p.add_argument("--gravity-scale", type=float, default=1.0)
    p.add_argument("--air-scale", type=float, default=1.0)
    p.add_argument("--rail-csv", type=Path, default=None, help="optional rail edge CSV")
    args = p.parse_args(argv)

    cfg = load_config()
    nodes = load_global_nodes(cfg)
    n_nodes = len(nodes)

    # Gravity matrices are independent of the airport cut-off: build them once.
    matrices = build_all_country_matrices(nodes, cfg, save=False)

    extra = None
    if args.rail_csv:
        from wwsim.rail import rail_mode_from_csv
        extra = [rail_mode_from_csv(args.rail_csv, nodes)]

    def _combine(air_edges: pd.DataFrame, stem: str) -> None:
        combined = build_combined_network(
            n_nodes, matrices, air_edges, cfg,
            gravity_scale=args.gravity_scale, air_scale=args.air_scale, extra_modes=extra,
        )
        save_combined_network(combined, cfg, stem=stem)
        logger.info("combined[%s]: %dx%d, %d nonzeros", stem, n_nodes, n_nodes, combined.nnz)

    if args.top_n is None:
        # Default: use the canonical air network produced by step 06.
        air_edges = pd.read_parquet(cfg.networks_dir / "global_air_network.parquet")
        _combine(air_edges, "global_combined_network")
    else:
        # Rebuild the air layer per cut-off from the step-03/04 interim artifacts.
        airport_edges = pd.read_parquet(cfg.interim_dir / "airport_edges.parquet")
        assignment = pd.read_parquet(cfg.interim_dir / "airport_node_assignment.parquet")
        multiple = len(args.top_n) > 1
        for n in args.top_n:
            cfg.air.top_n_airports = None if n == 0 else n
            label = "all" if cfg.air.top_n_airports is None else cfg.air.top_n_airports
            air_edges = build_air_network(airport_edges, assignment, cfg)
            stem = f"global_combined_network_top{label}" if multiple else "global_combined_network"
            _combine(air_edges, stem)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
