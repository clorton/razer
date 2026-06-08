"""OAG adapter: drop a licensed origin-destination passenger table into the pipeline.

The open OpenFlights + seat-capacity proxy (:mod:`wwsim.flights`) is a stand-in for OAG
because OAG passenger O-D volumes are commercial and have no free equivalent. When you do
have an OAG (or IATA AirportIS, Sabre, etc.) export, this adapter converts it into the
**same airport-edge schema** that :func:`wwsim.flights.build_airport_edges` produces, so
:func:`wwsim.air_network.build_air_network` consumes it with **zero downstream changes**.

Expected OAG CSV schema (one row per directed airport pair, e.g. annual totals):

    ``origin, destination, passengers``

where ``origin``/``destination`` are IATA codes. (Extra columns are ignored.) Country and
coordinates are joined from the airport master, exactly as the open path does.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .logging import logger


def load_oag_airport_edges(
    path: Path | str,
    airports: pd.DataFrame,
    origin_col: str = "origin",
    dest_col: str = "destination",
    pax_col: str = "passengers",
) -> pd.DataFrame:
    """Load a licensed OAG-style O-D table into the airport-edge schema.

    Args:
        path: Path to the OAG CSV.
        airports: Airport master (from :func:`wwsim.flights.load_airports`), providing
            ``iso3, lat, lon`` per IATA. Index must be IATA, or contain an ``iata`` column.
        origin_col: Column holding the origin IATA code.
        dest_col: Column holding the destination IATA code.
        pax_col: Column holding passenger volume (becomes ``weight``).

    Returns:
        DataFrame matching :func:`wwsim.flights.build_airport_edges`:
        ``[src_iata, dst_iata, src_iso3, dst_iso3, weight, n_carriers, src_lat, src_lon,
        dst_lat, dst_lon]`` (``n_carriers`` is set to 1; OAG O-D rows are pre-aggregated).

    Raises:
        KeyError: If required columns are missing from the OAG file.
    """
    df = pd.read_csv(path)
    for col in (origin_col, dest_col, pax_col):
        if col not in df.columns:
            raise KeyError(f"OAG file missing required column {col!r}")

    if airports.index.name != "iata" and "iata" in airports.columns:
        airports = airports.set_index("iata")

    edges = df.rename(
        columns={origin_col: "src_iata", dest_col: "dst_iata", pax_col: "weight"}
    )[["src_iata", "dst_iata", "weight"]].copy()
    edges = edges.groupby(["src_iata", "dst_iata"], as_index=False)["weight"].sum()
    edges["n_carriers"] = 1

    src = airports[["iso3", "lat", "lon"]].rename(
        columns={"iso3": "src_iso3", "lat": "src_lat", "lon": "src_lon"}
    )
    dst = airports[["iso3", "lat", "lon"]].rename(
        columns={"iso3": "dst_iso3", "lat": "dst_lat", "lon": "dst_lon"}
    )
    edges = edges.join(src, on="src_iata").join(dst, on="dst_iata")
    edges = edges.dropna(subset=["src_iso3", "dst_iso3"])
    logger.info("oag: loaded %d O-D airport edges from %s", len(edges), path)
    return edges
