"""Download open flight data and build weighted directed airport-to-airport edges.

Fetches OurAirports + OpenFlights, builds airport master and route edges with the
passenger-volume proxy, and writes intermediate artifacts for the air-network step.

Example:
    python scripts/03_ingest_flights.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim import flights  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402


def main() -> int:
    cfg = load_config()
    cfg.ensure_dirs()

    flights.download_flight_data(cfg)
    airports = flights.load_airports(cfg)
    routes = flights.load_routes(cfg)
    seats = flights.load_seat_table(cfg)
    edges = flights.build_airport_edges(routes, airports, seats, cfg)

    # Persist for downstream steps.
    airports.reset_index().to_parquet(cfg.interim_dir / "airports_master.parquet", index=False)
    edges.to_parquet(cfg.interim_dir / "airport_edges.parquet", index=False)

    # Summary.
    vol = flights.airport_volume(edges)
    cross = edges[edges["src_iso3"] != edges["dst_iso3"]]
    logger.info("airports=%d routes=%d edges=%d (cross-border=%d)",
                len(airports), len(routes), len(edges), len(cross))
    logger.info("airports missing ISO3: %d", int(airports["iso3"].isna().sum()))
    top = vol.head(10)
    logger.info("Top airports by passenger-proxy volume:\n%s", top.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
