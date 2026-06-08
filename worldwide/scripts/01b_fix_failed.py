"""Repair countries that failed acquisition, using the rasterio zonal-sum fallback.

The main acquisition (``01_acquire_countries.py``) uses RasterToolkit, which rejects
WorldPop rasters delivered on a global canvas (tie point x0 = -180) -- mostly antimeridian
crossers (RUS, FJI, KIR, TUV) plus a packaging quirk (ERI). This script re-runs only the
failed countries through :meth:`wwsim.acquire.Acquirer.acquire_country_rasterio` and updates
the manifest in place.

Example:
    python scripts/01b_fix_failed.py            # repair all failures in the manifest
    python scripts/01b_fix_failed.py RUS FJI    # repair specific countries
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.acquire import Acquirer  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402

MANIFEST_COLUMNS = ["iso3", "status", "source", "level", "n_units", "population", "gpkg", "error"]


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    cfg = load_config()
    manifest_path = cfg.nodes_dir / "acquisition_manifest.csv"
    df = pd.read_csv(manifest_path)
    manifest = df.set_index("iso3").to_dict("index")

    targets = [s.upper() for s in argv] if argv else [k for k, v in manifest.items() if v["status"] == "failed"]
    if not targets:
        logger.info("No failed countries to repair.")
        return 0

    logger.info("Repairing %d countries with rasterio fallback: %s", len(targets), targets)
    acq = Acquirer(cfg)
    for iso in targets:
        res = acq.acquire_country_rasterio(iso, force=True)
        manifest[iso] = {
            "status": res.status, "source": res.source, "level": res.level,
            "n_units": res.n_units, "population": res.population,
            "gpkg": str(res.gpkg) if res.gpkg else None, "error": res.error,
        }
        out = pd.DataFrame.from_dict(manifest, orient="index").reset_index(names="iso3")
        out.reindex(columns=MANIFEST_COLUMNS).to_csv(manifest_path, index=False)

    fixed = [t for t in targets if manifest[t]["status"] == "ok"]
    still = [t for t in targets if manifest[t]["status"] != "ok"]
    logger.info("Repair done: fixed=%s still_failed=%s", fixed, still)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
