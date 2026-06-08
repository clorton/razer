"""Tests for the pluggable rail and OAG adapters.

Failure here means a future rail dataset or a licensed OAG O-D table cannot be ingested
into the network without code changes -- the whole point of the adapters.
"""

from __future__ import annotations

import pandas as pd

from wwsim.oag import load_oag_airport_edges
from wwsim.rail import load_rail_edges_csv


def test_rail_direct_schema_passes_through(tmp_path):
    """Given a node-id rail CSV, when loaded, then the edges are returned unchanged."""
    csv = tmp_path / "rail.csv"
    pd.DataFrame(
        {"src_global_nodeid": [0, 1], "dst_global_nodeid": [2, 3], "weight": [5.0, 7.0]}
    ).to_csv(csv, index=False)
    edges = load_rail_edges_csv(csv)
    assert len(edges) == 2
    assert set(edges.columns) == {"src_global_nodeid", "dst_global_nodeid", "weight"}


def test_rail_coordinate_schema_snaps_to_nodes(tmp_path, global_nodes):
    """Given a coordinate rail CSV, when loaded, then endpoints snap to containing nodes.

    A cross-border rail link (AAA node 0 -> BBB node 2) must map to those global node ids.
    """
    csv = tmp_path / "rail_coords.csv"
    pd.DataFrame(
        {"src_lon": [0.5], "src_lat": [0.5], "dst_lon": [10.5], "dst_lat": [0.5], "weight": [9.0]}
    ).to_csv(csv, index=False)
    edges = load_rail_edges_csv(csv, global_nodes)
    assert len(edges) == 1
    assert edges.iloc[0]["src_global_nodeid"] == 0
    assert edges.iloc[0]["dst_global_nodeid"] == 2


def test_oag_loads_into_airport_edge_schema(tmp_path):
    """Given an OAG O-D CSV, when loaded, then it matches the airport-edge schema with ISO3.

    This is what lets a licensed OAG feed replace the open proxy with zero downstream change.
    """
    oag = tmp_path / "oag.csv"
    pd.DataFrame(
        {"origin": ["JFK", "JFK"], "destination": ["LHR", "LHR"], "passengers": [3000, 1000]}
    ).to_csv(oag, index=False)
    airports = pd.DataFrame(
        {"iata": ["JFK", "LHR"], "iso3": ["USA", "GBR"], "lat": [40.6, 51.5], "lon": [-73.8, -0.45]}
    ).set_index("iata")

    edges = load_oag_airport_edges(oag, airports)
    assert len(edges) == 1  # duplicate O-D rows summed
    row = edges.iloc[0]
    assert row["weight"] == 4000
    assert row["src_iso3"] == "USA" and row["dst_iso3"] == "GBR"
    assert {"src_lat", "src_lon", "dst_lat", "dst_lon", "n_carriers"}.issubset(edges.columns)
