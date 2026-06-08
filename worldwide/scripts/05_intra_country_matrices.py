"""Build per-country intra-country gravity migration matrices over admin-2 nodes.

Writes ``output/networks/<ISO3>_gravity.npz`` and ``<ISO3>_gravity.parquet`` per country.

Example:
    python scripts/05_intra_country_matrices.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.config import load_config  # noqa: E402
from wwsim.gravity import build_all_country_matrices  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def main() -> int:
    cfg = load_config()
    cfg.ensure_dirs()
    nodes = load_global_nodes(cfg)
    matrices = build_all_country_matrices(nodes, cfg, save=True)
    total_edges = sum(int((m > 0).sum()) for _, m in matrices.values())
    logger.info("gravity: %d countries, %d intra-country edges total", len(matrices), total_edges)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
