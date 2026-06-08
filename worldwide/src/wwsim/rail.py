"""Cross-border passenger rail: a documented, pluggable mode (not built by default).

Rail was a "nice to have". Rather than build a full rail ingester now, this module makes
rail a first-class :class:`wwsim.network.ModeNetwork` you can drop in later, and documents
the candidate open data sources and how they would map onto admin-2 nodes.

**Why air first, rail later.** For *worldwide inter-nation* spread, scheduled commercial
aviation dominates long-range mixing; cross-border rail matters in a few dense corridors
(continental Europe, and to a lesser degree parts of Asia and North America). So rail is an
additive refinement on specific corridors, not a global backbone -- exactly the kind of
thing the pluggable combiner is for.

**Candidate open sources** (all usable; effort/coverage trade-offs noted):

- **OpenStreetMap / OpenRailwayMap** (ODbL): global railway lines + stations. Cross-border
  links are derivable by intersecting ``railway=rail`` ways with country borders and
  snapping stations to admin-2. Highest coverage, heaviest processing (PBF extracts, e.g.
  via Geofabrik). Best long-term source.
- **Eurostat** ``rail_pa_*`` (e.g. ``rail_pa_intgonq``): EU international rail passenger
  flows by country pair -- ready-made cross-border *volumes* for Europe, country-level
  (would be distributed to admin-2 by population or by station location).
- **UIC Railway Statistics**: country-level totals; calibration only.
- **National open GTFS feeds** (DB, SNCF, Trenitalia, Amtrak, ...): scheduled services with
  station coordinates; international routes give explicit cross-border edges. Patchwork but
  precise where available.

**Mapping onto our node table.** Whatever the source, the deliverable is a directed edge
list ``[src_global_nodeid, dst_global_nodeid, weight]`` over admin-2 nodes -- produced by
snapping stations to the containing admin-2 (same point-in-polygon logic as airports, see
:mod:`wwsim.airports_assign`) and weighting by passenger volume (or service frequency as a
proxy). :func:`load_rail_edges_csv` accepts that, plus a convenience coordinate schema.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd

from .logging import logger
from .network import ModeNetwork


def load_rail_edges_csv(path: Path | str, global_nodes: gpd.GeoDataFrame | None = None) -> pd.DataFrame:
    """Load a user-provided rail edge table into the global node-id edge schema.

    Two input schemas are accepted:

    1. **Direct**: columns ``src_global_nodeid, dst_global_nodeid, weight`` -- returned as-is.
    2. **Coordinates**: columns ``src_lon, src_lat, dst_lon, dst_lat, weight`` -- each
       endpoint is mapped to the containing admin-2 node by point-in-polygon (requires
       ``global_nodes``).

    Args:
        path: Path to the rail edge CSV.
        global_nodes: Global node table (needed only for the coordinate schema).

    Returns:
        DataFrame ``[src_global_nodeid, dst_global_nodeid, weight]``.

    Raises:
        ValueError: If the CSV has neither recognized schema, or coordinates are given
            without ``global_nodes``.
    """
    df = pd.read_csv(path)
    direct = {"src_global_nodeid", "dst_global_nodeid", "weight"}
    coords = {"src_lon", "src_lat", "dst_lon", "dst_lat", "weight"}

    if direct.issubset(df.columns):
        return df[list(direct)].copy()

    if coords.issubset(df.columns):
        if global_nodes is None:
            raise ValueError("coordinate-schema rail edges require global_nodes for snapping")
        src = _snap_points(df["src_lon"], df["src_lat"], global_nodes)
        dst = _snap_points(df["dst_lon"], df["dst_lat"], global_nodes)
        out = pd.DataFrame(
            {"src_global_nodeid": src, "dst_global_nodeid": dst, "weight": df["weight"].to_numpy()}
        ).dropna()
        out["src_global_nodeid"] = out["src_global_nodeid"].astype("int64")
        out["dst_global_nodeid"] = out["dst_global_nodeid"].astype("int64")
        return out

    raise ValueError(
        "rail CSV must have either {src,dst}_global_nodeid+weight or {src,dst}_{lon,lat}+weight"
    )


def _snap_points(lon: pd.Series, lat: pd.Series, global_nodes: gpd.GeoDataFrame) -> pd.Series:
    """Map (lon, lat) points to the containing admin-2 global node id (NaN if none).

    Order and length match the input; points inside several overlapping polygons keep the
    first match.
    """
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon, lat), crs="EPSG:4326").reset_index(drop=True)
    joined = gpd.sjoin(
        pts, global_nodes[["global_nodeid", "geometry"]], predicate="within", how="left"
    )
    # sjoin keeps the left index; dedupe per original point, then restore input order.
    joined = joined[~joined.index.duplicated(keep="first")].sort_index()
    return joined["global_nodeid"].reset_index(drop=True)


def rail_mode_from_csv(
    path: Path | str, global_nodes: gpd.GeoDataFrame | None = None, scale: float = 1.0
) -> ModeNetwork:
    """Build a :class:`wwsim.network.ModeNetwork` for rail from a CSV.

    Args:
        path: Path to the rail edge CSV (see :func:`load_rail_edges_csv`).
        global_nodes: Global node table (for the coordinate schema).
        scale: Scale factor relative to other modes.

    Returns:
        A ``ModeNetwork`` named ``"rail"`` ready to pass to
        :func:`wwsim.network.combine_modes`.
    """
    edges = load_rail_edges_csv(path, global_nodes)
    logger.info("rail: loaded %d rail edges from %s", len(edges), path)
    return ModeNetwork("rail", edges, scale)
