"""Ingest open airport and air-route data, and build airport-to-airport edges.

This substitutes for commercial OAG passenger origin-destination data, which has no free
equivalent. We combine two open sources:

- **OurAirports** (``airports.csv``, public domain): the airport master -- coordinates,
  ISO country, IATA/ICAO codes. Clean ISO-2 country field for the cross-border test.
- **OpenFlights** (``routes.dat``, ODbL, ~2014 vintage): carrier-level origin-destination
  routes, with an ``equipment`` list of aircraft IATA type codes per route.

**Passenger-volume proxy.** OpenFlights has no seats or frequency. We weight each
carrier-route by a *seat-capacity* estimate: map each route's aircraft type code(s) to
typical seat counts (:data:`AIRCRAFT_SEATS`) and use their mean as that route's capacity.
Summing capacities across all carrier-routes on a directed airport pair gives a relative
"seats offered" volume -- a defensible stand-in for passenger volume. Alternatively the
``"routes"`` proxy just counts distinct carrier-routes (route multiplicity). A licensed OAG
O-D table can replace these weights via :mod:`wwsim.oag` without touching downstream code.

All edges here are **airport-level**; cross-border filtering, top-N selection, and
aggregation to admin-2 nodes happen in :mod:`wwsim.air_network`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pycountry
import requests

from .config import Config
from .logging import logger

OURAIRPORTS_URL = "https://davidmegginson.github.io/ourairports-data/airports.csv"
OPENFLIGHTS_BASE = "https://raw.githubusercontent.com/jpatokal/openflights/master/data/"
OPENFLIGHTS_FILES = {
    "routes": "routes.dat",
    "airports": "airports.dat",
    "planes": "planes.dat",
}

# OpenFlights routes.dat has no header; these are the documented column names, in order.
ROUTES_COLUMNS = [
    "airline", "airline_id", "src_iata", "src_id",
    "dst_iata", "dst_id", "codeshare", "stops", "equipment",
]

# Typical 1-2 class seat counts keyed by IATA aircraft type code (the codes that appear in
# OpenFlights routes.dat `equipment`). These are representative averages curated from
# manufacturer/airline typical-seating figures -- approximate by design, not authoritative.
# Unknown codes fall back to DEFAULT_SEATS (and are logged as a coverage gap).
AIRCRAFT_SEATS: dict[str, int] = {
    # --- Airbus narrowbody ---
    "318": 110, "319": 124, "31X": 124, "320": 150, "321": 185, "32A": 150, "32B": 185,
    "32C": 150, "32N": 165, "32Q": 195, "32S": 160, "323": 150,
    # --- Boeing 737 family ---
    "731": 110, "732": 120, "733": 128, "734": 147, "735": 110, "736": 110, "737": 130,
    "738": 162, "739": 178, "73C": 140, "73E": 140, "73G": 126, "73H": 162, "73J": 162,
    "73W": 126, "7M7": 162, "7M8": 178, "7M9": 193,
    # --- Boeing 747 ---
    "741": 366, "742": 366, "743": 400, "744": 416, "747": 416, "748": 410, "74E": 416,
    "74H": 416, "74L": 416, "74M": 366,
    # --- Boeing 757 / 767 ---
    "752": 200, "753": 243, "757": 200, "762": 245, "763": 269, "764": 245, "767": 250,
    "76W": 269,
    # --- Boeing 777 / 787 ---
    "772": 305, "773": 368, "77L": 317, "77W": 350, "777": 368, "778": 350, "779": 410,
    "787": 242, "788": 242, "789": 290, "78J": 330, "781": 330,
    # --- Airbus widebody ---
    "306": 266, "310": 220, "313": 220, "330": 290, "332": 253, "333": 295, "338": 287,
    "339": 300, "33X": 290, "340": 290, "342": 261, "343": 295, "345": 313, "346": 380,
    "350": 315, "351": 366, "358": 315, "359": 315, "35K": 366, "380": 525, "388": 525,
    # --- Regional jets (Bombardier/CRJ, Embraer, Fokker, etc.) ---
    "CR1": 50, "CR2": 50, "CR7": 70, "CR9": 90, "CRA": 90, "CRJ": 50, "CRK": 90,
    "E70": 70, "E75": 76, "E7W": 76, "E90": 96, "E95": 100, "EM2": 30, "EMB": 30,
    "ER3": 37, "ER4": 44, "ERD": 37, "ERJ": 50, "E170": 70, "E190": 100,
    "F70": 79, "F100": 100, "100": 100, "146": 100, "AR1": 100, "AR8": 100, "RJ1": 100,
    "RJ85": 100, "RJ100": 100, "SU9": 98, "SSJ": 98, "YK2": 120, "YK4": 120,
    # --- Turboprops / small ---
    "AT4": 48, "AT5": 48, "AT7": 70, "ATR": 70, "ATP": 64, "DH1": 37, "DH2": 37,
    "DH3": 50, "DH4": 78, "DH8": 50, "DHT": 19, "SF3": 34, "SFB": 34, "SWM": 19,
    "J31": 19, "J32": 19, "J41": 29, "BEH": 19, "BE1": 19, "BET": 9, "D38": 37,
    "EM2_": 30, "L4T": 19, "S20": 19, "SH6": 36, "BNI": 9, "CNA": 9, "CNC": 9, "PL2": 9,
    # --- MD / DC / older ---
    "M11": 285, "M1F": 285, "M80": 140, "M81": 140, "M82": 150, "M83": 150, "M87": 130,
    "M88": 150, "M90": 160, "D10": 270, "D11": 285, "D1C": 270, "D9S": 110, "D93": 120,
    "D94": 150, "D95": 120, "DC9": 110, "DC10": 270, "143": 100, "AB6": 260, "ABF": 260,
    "L10": 300, "IL9": 180, "I93": 180, "T20": 19, "T134": 80, "T154": 160,
}
DEFAULT_SEATS = 120  # fallback for aircraft type codes absent from AIRCRAFT_SEATS


def download_flight_data(cfg: Config, force: bool = False) -> dict[str, Path]:
    """Download the OurAirports master and OpenFlights route/airport/plane files.

    Args:
        cfg: Project configuration (provides ``flights_dir``).
        force: Re-download even if the file already exists.

    Returns:
        Mapping of logical name -> local path, e.g. ``{"airports_master": ..., "routes": ...}``.

    Raises:
        requests.HTTPError: If any download returns a non-200 status.
    """
    cfg.flights_dir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}

    targets = [("airports_master", OURAIRPORTS_URL, "ourairports.csv")]
    targets += [(name, OPENFLIGHTS_BASE + fn, fn) for name, fn in OPENFLIGHTS_FILES.items()]

    for name, url, local in targets:
        dest = cfg.flights_dir / local
        if dest.exists() and not force:
            logger.info("flights: %s already present", local)
            out[name] = dest
            continue
        logger.info("flights: downloading %s", url)
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        out[name] = dest

    # Also materialize the curated seat table for transparency/editing.
    seat_csv = cfg.flights_dir / "aircraft_seats.csv"
    if not seat_csv.exists() or force:
        pd.DataFrame(
            sorted(AIRCRAFT_SEATS.items()), columns=["type_code", "seats"]
        ).to_csv(seat_csv, index=False)
    out["seats"] = seat_csv
    return out


def alpha2_to_alpha3(alpha2: str) -> str | None:
    """Convert an ISO 3166-1 alpha-2 code to alpha-3.

    Args:
        alpha2: Two-letter country code (as used by OurAirports ``iso_country``).

    Returns:
        The alpha-3 code, or ``None`` if unknown (e.g. ``"XK"`` Kosovo, not a UN member).
    """
    if not isinstance(alpha2, str) or len(alpha2) != 2:
        return None
    rec = pycountry.countries.get(alpha_2=alpha2.upper())
    return rec.alpha_3 if rec else None


def load_airports(cfg: Config) -> pd.DataFrame:
    """Load the airport master keyed by IATA code, with ISO3 country and coordinates.

    Built from OurAirports; rows without an IATA code or coordinates are dropped. When a
    single IATA code maps to several airport records, the one with the most "serious" type
    (large > medium > small > other) and scheduled service is kept.

    Args:
        cfg: Project configuration.

    Returns:
        DataFrame indexed by ``iata`` with columns
        ``[icao, name, iso2, iso3, lat, lon, type, scheduled]``.
    """
    path = cfg.flights_dir / "ourairports.csv"
    df = pd.read_csv(path, low_memory=False)

    df = df[df["iata_code"].notna() & (df["iata_code"].str.len() == 3)].copy()
    df = df.dropna(subset=["latitude_deg", "longitude_deg"])
    df["iso3"] = df["iso_country"].map(alpha2_to_alpha3)

    type_rank = {"large_airport": 0, "medium_airport": 1, "small_airport": 2}
    df["_rank"] = df["type"].map(type_rank).fillna(3)
    df["_sched"] = (df.get("scheduled_service", "no") == "yes").astype(int)
    # Prefer scheduled, then more serious type.
    df = df.sort_values(["iata_code", "_sched", "_rank"], ascending=[True, False, True])
    df = df.drop_duplicates("iata_code", keep="first")

    out = pd.DataFrame({
        "iata": df["iata_code"].values,
        "icao": df["icao_code"].values,
        "name": df["name"].values,
        "iso2": df["iso_country"].values,
        "iso3": df["iso3"].values,
        "lat": df["latitude_deg"].astype(float).values,
        "lon": df["longitude_deg"].astype(float).values,
        "type": df["type"].values,
        "scheduled": df["_sched"].astype(bool).values,
    }).set_index("iata")
    logger.info("flights: loaded %d airports (with IATA + coords)", len(out))
    return out


def load_routes(cfg: Config) -> pd.DataFrame:
    """Load OpenFlights routes, parsing the ``\\N`` null sentinel and equipment list.

    Args:
        cfg: Project configuration.

    Returns:
        DataFrame with columns ``[airline, src_iata, dst_iata, codeshare, stops, equipment]``
        where ``codeshare`` is bool, ``stops`` is int, and ``equipment`` is a list of
        aircraft type codes.
    """
    path = cfg.flights_dir / "routes.dat"
    df = pd.read_csv(path, header=None, names=ROUTES_COLUMNS, na_values=["\\N"], keep_default_na=True)

    df["codeshare"] = df["codeshare"].fillna("").astype(str).str.upper().eq("Y")
    df["stops"] = pd.to_numeric(df["stops"], errors="coerce").fillna(0).astype(int)
    df["equipment"] = (
        df["equipment"].fillna("").astype(str).str.split().apply(lambda xs: [x for x in xs if x])
    )
    df = df.dropna(subset=["src_iata", "dst_iata"])
    df = df[(df["src_iata"].str.len() == 3) & (df["dst_iata"].str.len() == 3)]
    logger.info("flights: loaded %d routes", len(df))
    return df[["airline", "src_iata", "dst_iata", "codeshare", "stops", "equipment"]]


def load_seat_table(cfg: Config) -> dict[str, int]:
    """Load the aircraft-type -> seats lookup (curated default, overridable via CSV).

    Args:
        cfg: Project configuration.

    Returns:
        Mapping of aircraft IATA type code to typical seat count.
    """
    seat_csv = cfg.flights_dir / "aircraft_seats.csv"
    if seat_csv.exists():
        df = pd.read_csv(seat_csv)
        return dict(zip(df["type_code"].astype(str), df["seats"].astype(int)))
    return dict(AIRCRAFT_SEATS)


def _route_capacity(equipment: list[str], seats: dict[str, int]) -> tuple[float, int]:
    """Representative seat capacity for one carrier-route and its type-code coverage.

    Args:
        equipment: List of aircraft IATA type codes operating the route.
        seats: Type-code -> seats lookup.

    Returns:
        Tuple of (mean seats across listed types, number of listed types resolved by the
        lookup). When no type is given, returns ``(DEFAULT_SEATS, 0)``.
    """
    if not equipment:
        return float(DEFAULT_SEATS), 0
    vals = [seats.get(code, DEFAULT_SEATS) for code in equipment]
    resolved = sum(1 for code in equipment if code in seats)
    return float(np.mean(vals)), resolved


def build_airport_edges(
    routes: pd.DataFrame,
    airports: pd.DataFrame,
    seats: dict[str, int],
    cfg: Config,
) -> pd.DataFrame:
    """Aggregate carrier-routes into weighted directed airport-to-airport edges.

    Applies the configured filters (codeshare/direct), joins airports to obtain ISO3
    countries and coordinates, weights each carrier-route by the passenger proxy, and sums
    over carriers per directed airport pair.

    Args:
        routes: Output of :func:`load_routes`.
        airports: Output of :func:`load_airports`.
        seats: Output of :func:`load_seat_table`.
        cfg: Project configuration (uses ``cfg.air``).

    Returns:
        DataFrame with columns
        ``[src_iata, dst_iata, src_iso3, dst_iso3, weight, n_carriers, src_lat, src_lon,
        dst_lat, dst_lon]``. ``weight`` is summed seats (``passenger_proxy="seats"``) or
        carrier-route count (``"routes"``).
    """
    ap = cfg.air
    df = routes
    if ap.exclude_codeshares:
        df = df[~df["codeshare"]]
    if ap.direct_only:
        df = df[df["stops"] == 0]

    # Keep only routes whose endpoints exist in the airport master.
    valid = airports.index
    df = df[df["src_iata"].isin(valid) & df["dst_iata"].isin(valid)].copy()

    if ap.passenger_proxy == "seats":
        caps = df["equipment"].apply(lambda e: _route_capacity(e, seats))
        df["weight"] = [c for c, _ in caps]
        resolved = sum(r for _, r in caps)
        total_codes = int(df["equipment"].apply(len).sum())
        if total_codes:
            logger.info(
                "flights: seat-table covers %.1f%% of route equipment codes",
                100.0 * resolved / total_codes,
            )
    elif ap.passenger_proxy == "routes":
        df["weight"] = 1.0
    else:
        raise ValueError(f"unknown passenger_proxy {ap.passenger_proxy!r}")

    grouped = (
        df.groupby(["src_iata", "dst_iata"], as_index=False)
        .agg(weight=("weight", "sum"), n_carriers=("airline", "count"))
    )

    # Attach country + coordinates for both endpoints.
    src = airports[["iso3", "lat", "lon"]].rename(
        columns={"iso3": "src_iso3", "lat": "src_lat", "lon": "src_lon"}
    )
    dst = airports[["iso3", "lat", "lon"]].rename(
        columns={"iso3": "dst_iso3", "lat": "dst_lat", "lon": "dst_lon"}
    )
    grouped = grouped.join(src, on="src_iata").join(dst, on="dst_iata")
    logger.info("flights: built %d directed airport edges", len(grouped))
    return grouped


def airport_volume(edges: pd.DataFrame) -> pd.Series:
    """Total passenger-proxy throughput per airport (in + out), for top-N ranking.

    Args:
        edges: Output of :func:`build_airport_edges`.

    Returns:
        Series indexed by IATA code, summed inbound + outbound weight, descending.
    """
    out = edges.groupby("src_iata")["weight"].sum()
    inb = edges.groupby("dst_iata")["weight"].sum()
    vol = out.add(inb, fill_value=0.0).sort_values(ascending=False)
    vol.index.name = "iata"
    return vol
