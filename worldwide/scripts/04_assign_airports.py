"""Assign airports to admin-2 nodes (point-in-polygon) and persist the mapping.

Reads the global node table and the airport master, writes
``data/interim/airport_node_assignment.parquet``.

Example:
    python scripts/04_assign_airports.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.airports_assign import assign_airports_to_nodes  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def main() -> int:
    cfg = load_config()
    nodes = load_global_nodes(cfg)
    airports = pd.read_parquet(cfg.interim_dir / "airports_master.parquet")
    assignment = assign_airports_to_nodes(airports, nodes)
    out = cfg.interim_dir / "airport_node_assignment.parquet"
    assignment.to_parquet(out, index=False)
    logger.info("assignment: wrote %s (%d airports)", out, len(assignment))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
