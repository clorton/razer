"""Run the full post-acquisition pipeline (steps 02-08) in order.

Assumes steps 01 (country acquisition) and 03 (flight ingestion) inputs exist, or runs them
if requested. Each step is idempotent and reads/writes the on-disk layout from
:mod:`wwsim.config`.

Example:
    python scripts/run_all.py             # 02 -> 08 (nodes, assign, matrices, networks, plots)
    python scripts/run_all.py --with-acquire --with-flights   # also run 01 and 03 first
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS.parents[0] / "src"))

from wwsim.logging import logger  # noqa: E402


def _run(script_name: str, argv: list[str] | None = None) -> None:
    """Import a sibling script module and call its main()."""
    path = SCRIPTS / script_name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    logger.info("=== running %s ===", script_name)
    rc = mod.main(argv) if argv is not None else mod.main()
    if rc:
        raise SystemExit(rc)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--with-acquire", action="store_true", help="run 01 (acquire countries) first")
    p.add_argument("--with-flights", action="store_true", help="run 03 (ingest flights) first")
    args = p.parse_args()

    if args.with_acquire:
        _run("01_acquire_countries.py", [])
        _run("01b_fix_failed.py", [])
    _run("02_build_nodes.py")
    if args.with_flights:
        _run("03_ingest_flights.py")
    _run("04_assign_airports.py")
    _run("05_intra_country_matrices.py")
    _run("06_global_air_network.py")
    _run("07_combine_network.py", [])
    _run("08_plots.py", [])
    logger.info("run_all: complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
