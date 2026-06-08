"""Fuse per-country admin-2 GeoPackages into the unified global node table.

Reads every ``data/countries/<ISO3>/<ISO3>_admin*.gpkg`` and writes
``output/nodes/global_admin2_nodes.{gpkg,parquet}``.

Example:
    python scripts/02_build_nodes.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import build_global_nodes, save_global_nodes  # noqa: E402


def main() -> int:
    cfg = load_config()
    cfg.ensure_dirs()
    nodes = build_global_nodes(cfg)
    save_global_nodes(nodes, cfg)
    logger.info(
        "nodes: %d countries, %d admin units, total pop=%s",
        nodes["iso3"].nunique(), len(nodes), f"{nodes['population'].sum():,.0f}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
