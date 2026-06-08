"""Acquire admin-2 boundaries + 1km WorldPop population for UN member states.

Drives :class:`wwsim.acquire.Acquirer` over every UN member (or a subset given on the
command line), writing one ``<ISO3>_admin2.gpkg`` per country and an incremental manifest
so the run is fully resumable -- re-running skips countries whose GeoPackage already exists.

Examples:
    # All 193 UN members, default source waterfall (unocha -> geoboundaries -> gadm):
    python scripts/01_acquire_countries.py

    # A subset, forcing geoBoundaries first (fast per-country downloads, no 1-2GB UNOCHA):
    python scripts/01_acquire_countries.py NGA KEN FRA --source-order geoboundaries gadm

    # Re-acquire a country from scratch:
    python scripts/01_acquire_countries.py COM --force
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Allow running as a plain script (python scripts/01_...) without installing on PYTHONPATH.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.acquire import Acquirer  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.countries import UN_MEMBERS  # noqa: E402
from wwsim.logging import logger  # noqa: E402

MANIFEST_COLUMNS = ["iso3", "status", "source", "level", "n_units", "population", "gpkg", "error"]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("isos", nargs="*", help="ISO3 codes to acquire (default: all 193 UN members)")
    p.add_argument("--source-order", nargs="+", default=None,
                   help="Override shape-source waterfall, e.g. --source-order geoboundaries gadm")
    p.add_argument("--force", action="store_true", help="Re-acquire even if a gpkg exists")
    p.add_argument("--limit", type=int, default=None, help="Process at most this many countries")
    p.add_argument("--config", type=Path, default=None, help="Optional YAML config overlay")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cfg = load_config(args.config)
    if args.source_order:
        cfg.shape_source_order = tuple(s.lower() for s in args.source_order)
    cfg.ensure_dirs()

    isos = [s.upper() for s in args.isos] if args.isos else list(UN_MEMBERS)
    if args.limit:
        isos = isos[: args.limit]

    manifest_path = cfg.nodes_dir / "acquisition_manifest.csv"
    # Resume from an existing manifest so re-runs append rather than overwrite.
    if manifest_path.exists() and not args.force:
        manifest = pd.read_csv(manifest_path).set_index("iso3").to_dict("index")
    else:
        manifest = {}

    logger.info("Acquiring %d countries; order=%s", len(isos), cfg.shape_source_order)
    acq = Acquirer(cfg)

    ok = skipped = failed = 0
    for i, iso in enumerate(isos, 1):
        logger.info("[%d/%d] %s", i, len(isos), iso)
        res = acq.acquire_country(iso, force=args.force)
        manifest[iso] = {
            "status": res.status, "source": res.source, "level": res.level,
            "n_units": res.n_units, "population": res.population,
            "gpkg": str(res.gpkg) if res.gpkg else None, "error": res.error,
        }
        ok += res.status == "ok"
        skipped += res.status == "skipped"
        failed += res.status == "failed"

        # Persist after every country so an interrupted run loses nothing.
        df = pd.DataFrame.from_dict(manifest, orient="index").reset_index(names="iso3")
        df.reindex(columns=MANIFEST_COLUMNS).to_csv(manifest_path, index=False)

    logger.info("DONE: ok=%d skipped=%d failed=%d -> %s", ok, skipped, failed, manifest_path)
    failures = [k for k, v in manifest.items() if v["status"] == "failed"]
    if failures:
        logger.warning("Failed: %s", ", ".join(sorted(failures)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
