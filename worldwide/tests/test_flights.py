"""Tests for flight data parsing and airport-edge construction.

Failure here means the air network is built from corrupted routes or a broken passenger
proxy, distorting inter-country coupling.
"""

from __future__ import annotations

import pandas as pd

from wwsim.config import Config
from wwsim.flights import (
    _route_capacity,
    airport_volume,
    alpha2_to_alpha3,
    build_airport_edges,
    load_routes,
)


def test_alpha2_to_alpha3_maps_and_handles_unknown():
    """Given ISO-2 codes, when converted, then valid ones map and unknowns return None."""
    assert alpha2_to_alpha3("US") == "USA"
    assert alpha2_to_alpha3("FR") == "FRA"
    assert alpha2_to_alpha3("ZZ") is None


def test_route_capacity_uses_seat_table_mean():
    """Given equipment codes, when capacity computed, then it is the mean of looked-up seats."""
    seats = {"320": 150, "738": 160}
    cap, resolved = _route_capacity(["320", "738"], seats)
    assert cap == 155.0
    assert resolved == 2
    # Unknown code falls back to default and is not counted as resolved.
    cap2, resolved2 = _route_capacity(["ZZZ"], seats)
    assert resolved2 == 0


def test_load_routes_parses_null_sentinel_and_equipment(tmp_path):
    """Given a routes.dat with \\N nulls, when loaded, then nulls/equipment parse correctly.

    OpenFlights uses the literal '\\N' for missing values; mis-parsing it would poison joins.
    """
    cfg = Config(project_root=tmp_path)
    cfg.ensure_dirs()
    (cfg.flights_dir / "routes.dat").write_text(
        "AA,24,JFK,3797,LHR,507,,0,320 321\n"
        "BB,\\N,XXX,\\N,YYY,\\N,Y,0,738\n"
    )
    routes = load_routes(cfg)
    assert len(routes) == 2
    assert routes.iloc[0]["equipment"] == ["320", "321"]
    assert routes.iloc[0]["codeshare"] is False or routes.iloc[0]["codeshare"] == False  # noqa: E712
    assert routes.iloc[1]["codeshare"] == True  # noqa: E712


def test_build_airport_edges_filters_and_weights():
    """Given routes + airports, when edges built, then codeshares drop and seats sum per pair.

    Validates the seats proxy and the codeshare/direct filters that prevent double counting.
    """
    cfg = Config()
    cfg.air.passenger_proxy = "seats"
    cfg.air.exclude_codeshares = True
    cfg.air.direct_only = True

    airports = pd.DataFrame(
        {"iata": ["JFK", "LHR"], "iso3": ["USA", "GBR"], "lat": [40.6, 51.5], "lon": [-73.8, -0.45]}
    ).set_index("iata")
    routes = pd.DataFrame(
        {
            "airline": ["AA", "BA", "VS"],
            "src_iata": ["JFK", "JFK", "JFK"],
            "dst_iata": ["LHR", "LHR", "LHR"],
            "codeshare": [False, False, True],  # the VS row is a codeshare -> dropped
            "stops": [0, 0, 0],
            "equipment": [["77W"], ["744"], ["320"]],
        }
    )
    seats = {"77W": 350, "744": 416, "320": 150}
    edges = build_airport_edges(routes, airports, seats, cfg)
    assert len(edges) == 1
    row = edges.iloc[0]
    assert row["src_iata"] == "JFK" and row["dst_iata"] == "LHR"
    assert row["n_carriers"] == 2  # codeshare excluded
    assert row["weight"] == 350 + 416  # summed seats of the two non-codeshare carriers
    assert row["src_iso3"] == "USA" and row["dst_iso3"] == "GBR"


def test_airport_volume_sums_in_and_out():
    """Given directed edges, when volume computed, then it sums inbound + outbound weight."""
    edges = pd.DataFrame(
        {"src_iata": ["A", "B"], "dst_iata": ["B", "A"], "weight": [10.0, 4.0]}
    )
    vol = airport_volume(edges)
    assert vol.loc["A"] == 14.0
    assert vol.loc["B"] == 14.0
