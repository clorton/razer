"""Fetch real COVID-19 data (Our World in Data) and plot it for comparison with the sim.

COVID-19 ran ~2020-2024 (it did not exist in 2000-2004), so this pulls the 2020-2024 OWID
cases/deaths series. Writes plots into ``output/covid/``; if a sim run exists
(``output/seir/totals_*.csv``) it also overlays the single-wave shape.

Example:
    python scripts/11_covid_reference.py
    python scripts/11_covid_reference.py --sim-totals output/seir/totals_demo_air.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim import covid_reference as cov  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402

# OWID location names matching the sim's sample ISO3 set (CHN/USA/GBR/ITA/BRA/ZAF/IND/AUS).
SAMPLE_COUNTRIES = ["China", "United States", "United Kingdom", "Italy",
                    "Brazil", "South Africa", "India", "Australia"]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--force", action="store_true", help="re-download the OWID CSV")
    p.add_argument("--sim-totals", type=Path, default=Path("output/seir/totals_demo_air.csv"),
                   help="a wwsim SEIR totals CSV for the shape overlay (optional)")
    args = p.parse_args(argv)

    cfg = load_config()
    out_dir = cfg.output_dir / "covid"
    out_dir.mkdir(parents=True, exist_ok=True)

    cov.download_owid(cfg, force=args.force)
    df = cov.load_owid(cfg)

    w = cov.location_series(df, "World")
    logger.info("covid: World total cases=%s, total deaths=%s, dates %s..%s",
                f"{w['total_cases'].max():,.0f}", f"{w['total_deaths'].max():,.0f}",
                w.index.min().date(), w.index.max().date())

    cov.plot_global(df, out_dir / "covid_global.png")
    cov.plot_countries(df, SAMPLE_COUNTRIES, out_dir / "covid_countries.png")

    if args.sim_totals.exists():
        sim = pd.read_csv(args.sim_totals)
        cov.plot_first_wave_vs_sim(df, sim, out_dir / "covid_vs_sim_firstwave.png")
    else:
        logger.info("covid: %s not found; skipping sim overlay", args.sim_totals)

    logger.info("covid: outputs -> %s", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
