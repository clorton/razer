"""Tests for airport -> admin-2 node assignment.

Failure here means airports attach to the wrong node (or none), corrupting which admin-2
units the air network connects across borders.
"""

from __future__ import annotations

import pandas as pd

from wwsim.airports_assign import assign_airports_to_nodes


def test_within_assignment(global_nodes):
    """Given an airport inside a node, when assigned, then it maps to that node via 'within'."""
    airports = pd.DataFrame(
        {"iata": ["AA1"], "lon": [0.5], "lat": [0.5], "iso3": ["AAA"]}
    )
    out = assign_airports_to_nodes(airports, global_nodes)
    row = out[out["iata"] == "AA1"].iloc[0]
    assert row["node_global_nodeid"] == 0
    assert row["node_iso3"] == "AAA"
    assert row["assign_method"] == "within"


def test_nearest_fallback_for_offshore_airport(global_nodes):
    """Given an airport outside all polygons but near one, when assigned, then 'nearest' snaps it.

    Coastal/offshore airports often miss containment; the nearest-in-country fallback must
    still attach them to their own country's node.
    """
    # lon 1.2 sits in the gap between AAA's two squares (1..2), within 0.5deg of node 0.
    airports = pd.DataFrame(
        {"iata": ["OFF"], "lon": [1.2], "lat": [0.5], "iso3": ["AAA"]}
    )
    out = assign_airports_to_nodes(airports, global_nodes)
    row = out[out["iata"] == "OFF"].iloc[0]
    assert row["assign_method"] == "nearest"
    assert row["node_iso3"] == "AAA"
    assert row["node_global_nodeid"] == 0  # node 0 (edge at lon 1) is nearest


def test_border_disambiguation_prefers_matching_country(global_nodes):
    """Given overlapping border polygons, when assigned, then the airport's own country wins.

    Different sources can make neighboring countries' polygons overlap; the airport should
    bind to the node whose ISO3 matches the airport, keeping the cross-border test correct.
    """
    import geopandas as gpd
    from shapely.geometry import box

    # Add a BBB node that overlaps AAA's node 0 region (a deliberate border overlap).
    overlap = gpd.GeoDataFrame(
        [[4, "BBB", "B-overlap", 50.0, 0.5, 0.5, box(0, 0, 1, 1)]],
        columns=["global_nodeid", "iso3", "adm2_name", "population", "lon", "lat", "geometry"],
        geometry="geometry", crs="EPSG:4326",
    )
    nodes = gpd.GeoDataFrame(pd.concat([global_nodes, overlap], ignore_index=True),
                             geometry="geometry", crs="EPSG:4326")
    airports = pd.DataFrame({"iata": ["AA1"], "lon": [0.5], "lat": [0.5], "iso3": ["AAA"]})
    out = assign_airports_to_nodes(airports, nodes)
    row = out[out["iata"] == "AA1"].iloc[0]
    assert row["node_iso3"] == "AAA"  # not the overlapping BBB node
    assert row["node_global_nodeid"] == 0
