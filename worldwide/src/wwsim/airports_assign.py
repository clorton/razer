"""Assign each airport to the admin-2 node that contains it.

The global air network links admin-2 *nodes* (not airports), so every airport must be
mapped to a node. We do this by point-in-polygon containment against the global node
table, with two refinements:

1. **Border disambiguation.** Different countries' admin-2 polygons (drawn from different
   sources) can overlap slightly at borders, so a point may fall "within" more than one
   node. When that happens we keep the node whose country matches the airport's own ISO3.
2. **Nearest fallback.** Airports just offshore, on reclaimed land, or in a country we
   could not acquire may not fall within any polygon. We snap them to the nearest node in
   the airport's own country (within a distance cap); failing that, to the nearest node
   globally. Unassigned airports are logged and dropped from the network.

The result keys each airport (IATA) to a ``node_global_nodeid`` and that node's ISO3,
which the air network uses for the cross-border test and admin-2 aggregation.
"""

from __future__ import annotations

import warnings

import geopandas as gpd
import pandas as pd

from .logging import logger

# Max snap distance for the nearest-node fallback, in degrees (~0.5 deg ~ 55 km at equator).
_MAX_SNAP_DEG = 0.5


def _points_gdf(airports: pd.DataFrame) -> gpd.GeoDataFrame:
    """Build a point GeoDataFrame from an airport table with lon/lat columns."""
    return gpd.GeoDataFrame(
        airports.copy(),
        geometry=gpd.points_from_xy(airports["lon"], airports["lat"]),
        crs="EPSG:4326",
    )


def assign_airports_to_nodes(
    airports: pd.DataFrame, global_nodes: gpd.GeoDataFrame
) -> pd.DataFrame:
    """Map airports to admin-2 nodes by containment, with same-country nearest fallback.

    Args:
        airports: Airport master with at least ``iata, lon, lat, iso3`` (the airport's own
            country from OurAirports). Typically the output of
            :func:`wwsim.flights.load_airports` (reset to a column index).
        global_nodes: Global node table with ``global_nodeid, iso3, geometry`` (EPSG:4326).

    Returns:
        DataFrame with columns ``[iata, airport_iso3, node_global_nodeid, node_iso3,
        assign_method]`` for every airport that could be assigned. ``assign_method`` is
        ``"within"`` or ``"nearest"``.
    """
    if "iata" not in airports.columns:
        airports = airports.reset_index()

    pts = _points_gdf(airports[["iata", "lon", "lat", "iso3"]].rename(columns={"iso3": "airport_iso3"}))
    nodes = global_nodes[["global_nodeid", "iso3", "geometry"]].rename(columns={"iso3": "node_iso3"})

    # --- 1) containment join ---
    within = gpd.sjoin(pts, nodes, predicate="within", how="left")
    within = within.rename(columns={"global_nodeid": "node_global_nodeid"})

    # Disambiguate multi-matches: prefer the node whose country matches the airport.
    within["_match"] = (within["airport_iso3"] == within["node_iso3"]).astype(int)
    within = (
        within.sort_values(["iata", "_match"], ascending=[True, False])
        .drop_duplicates("iata", keep="first")
    )

    assigned = within[within["node_global_nodeid"].notna()].copy()
    assigned["assign_method"] = "within"

    # --- 2) nearest fallback for the unmatched ---
    unmatched_iatas = set(pts["iata"]) - set(assigned["iata"])
    fallback_rows: list[dict] = []
    if unmatched_iatas:
        un = pts[pts["iata"].isin(unmatched_iatas)]
        nodes_proj = nodes  # already 4326; sjoin_nearest works in CRS units (degrees)
        # Nearest within the airport's own country first (keeps cross-border test correct).
        for airport_iso3, grp in un.groupby("airport_iso3"):
            country_nodes = nodes_proj[nodes_proj["node_iso3"] == airport_iso3]
            if len(country_nodes) == 0:
                continue
            with warnings.catch_warnings():
                # Nearest in a geographic CRS is approximate; fine for this short-range snap.
                warnings.simplefilter("ignore")
                near = gpd.sjoin_nearest(
                    grp, country_nodes, how="left", max_distance=_MAX_SNAP_DEG,
                    distance_col="_dist",
                ).rename(columns={"global_nodeid": "node_global_nodeid"})
            near = near.dropna(subset=["node_global_nodeid"]).drop_duplicates("iata")
            for _, r in near.iterrows():
                fallback_rows.append({
                    "iata": r["iata"], "airport_iso3": r["airport_iso3"],
                    "node_global_nodeid": r["node_global_nodeid"], "node_iso3": r["node_iso3"],
                    "assign_method": "nearest",
                })

    out = pd.concat(
        [
            assigned[["iata", "airport_iso3", "node_global_nodeid", "node_iso3", "assign_method"]],
            pd.DataFrame(fallback_rows),
        ],
        ignore_index=True,
    )
    out["node_global_nodeid"] = out["node_global_nodeid"].astype("int64")

    n_total = len(pts)
    n_within = int((out["assign_method"] == "within").sum())
    n_near = int((out["assign_method"] == "nearest").sum())
    logger.info(
        "airports_assign: %d/%d assigned (within=%d, nearest=%d, unassigned=%d)",
        len(out), n_total, n_within, n_near, n_total - len(out),
    )
    return out
