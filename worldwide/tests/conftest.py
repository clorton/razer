"""Shared pytest fixtures: a temp Config and a small synthetic global node table.

The synthetic world has two countries (``AAA``, ``BBB``), each with two square admin-2
units laid out along the equator, so containment, distance, and cross-border logic are all
easy to reason about by hand.
"""

from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim.config import Config  # noqa: E402


@pytest.fixture
def cfg(tmp_path) -> Config:
    """A Config rooted at a temp dir so tests never touch real data/output."""
    c = Config(project_root=tmp_path)
    c.ensure_dirs()
    return c


@pytest.fixture
def global_nodes() -> gpd.GeoDataFrame:
    """Four admin-2 nodes across two countries.

    Layout (each a 1x1 deg square; centroids in parentheses):
        AAA: node 0 lon[0,1] (0.5, 0.5), node 1 lon[2,3] (2.5, 0.5)
        BBB: node 2 lon[10,11] (10.5, 0.5), node 3 lon[12,13] (12.5, 0.5)
    """
    rows = [
        (0, "AAA", "A-west", 100.0, 0.5, 0.5, box(0, 0, 1, 1)),
        (1, "AAA", "A-east", 200.0, 2.5, 0.5, box(2, 0, 3, 1)),
        (2, "BBB", "B-west", 300.0, 10.5, 0.5, box(10, 0, 11, 1)),
        (3, "BBB", "B-east", 400.0, 12.5, 0.5, box(12, 0, 13, 1)),
    ]
    gdf = gpd.GeoDataFrame(
        rows,
        columns=["global_nodeid", "iso3", "adm2_name", "population", "lon", "lat", "geometry"],
        geometry="geometry",
        crs="EPSG:4326",
    )
    return gdf
