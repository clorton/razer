"""Tests for the unified node-table normalizer.

Failure here means the heterogeneous per-source GeoPackages are not harmonized, so the
global node table loses ids, names, or population.
"""

from __future__ import annotations

import geopandas as gpd
from shapely.geometry import box

from wwsim.nodes import fix_mojibake, normalize_country_nodes


def test_fix_mojibake_repairs_double_encoding():
    """Given a mojibake string, when fixed, then the accented form is restored."""
    assert fix_mojibake("MboudÃ©") == "Mboudé"  # "MboudÃ©" -> "Mboudé"
    assert fix_mojibake("Lagos") == "Lagos"  # untouched


def test_normalize_geoboundaries_schema(tmp_path):
    """Given a geoBoundaries-style gpkg, when normalized, then shapeID/shapeName map through.

    geoBoundaries gpkgs carry shapeName/shapeID and no iso3; the normalizer must inject iso3
    and key the unified schema off those columns.
    """
    gdf = gpd.GeoDataFrame(
        {
            "shapeName": ["North", "South"],
            "shapeID": ["X-1", "X-2"],
            "nodeid": [0, 1],
            "name": ["North", "South"],
            "population": [1000, None],  # a null population must become 0
            "geometry": [box(0, 0, 1, 1), box(1, 0, 2, 1)],
        },
        geometry="geometry", crs="EPSG:4326",
    )
    path = tmp_path / "XYZ_admin2.gpkg"
    gdf.to_file(path, driver="GPKG")

    out = normalize_country_nodes("XYZ", path)
    assert list(out["iso3"].unique()) == ["XYZ"]
    assert set(out["adm2_id"]) == {"X-1", "X-2"}
    assert out["population"].sum() == 1000.0  # null -> 0
    assert {"lon", "lat", "adm2_name", "geometry"}.issubset(out.columns)


def test_normalize_unocha_schema_uses_pcode_and_adm1(tmp_path):
    """Given a UNOCHA-style gpkg, when normalized, then adm2_pcode/adm1_name are picked up."""
    gdf = gpd.GeoDataFrame(
        {
            "adm0_name": ["Country", "Country"],
            "adm1_name": ["Prov1", "Prov1"],
            "adm2_name": ["D1", "D2"],
            "adm2_pcode": ["PC1", "PC2"],
            "nodeid": [0, 1],
            "name": ["D1", "D2"],
            "population": [10, 20],
            "geometry": [box(0, 0, 1, 1), box(1, 0, 2, 1)],
        },
        geometry="geometry", crs="EPSG:4326",
    )
    path = tmp_path / "ABC_admin2.gpkg"
    gdf.to_file(path, driver="GPKG")

    out = normalize_country_nodes("ABC", path)
    assert set(out["adm2_id"]) == {"PC1", "PC2"}
    assert list(out["adm1_name"].unique()) == ["Prov1"]
