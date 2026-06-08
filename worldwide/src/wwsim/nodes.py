"""Fuse per-country admin-2 GeoPackages into one global node table.

Each country's GeoPackage (from :mod:`wwsim.acquire`) carries source-specific columns
(UNOCHA ``adm2_pcode``, geoBoundaries ``shapeID``, GADM ``GID_2`` -- with admin-1 variants
when a country falls back a level). This module normalizes them to a single schema and
stacks them into one worldwide GeoDataFrame where every admin-2 unit gets a stable
``global_nodeid``. That node table is the backbone of every network: per-country gravity
matrices index into it, and the global air network connects its nodes across borders.

Unified node schema:
    ``global_nodeid, iso3, adm2_name, adm2_id, adm1_name, population, lon, lat, geometry``
"""

from __future__ import annotations

import warnings
from pathlib import Path

import geopandas as gpd
import pandas as pd

from .config import Config
from .logging import logger

# Candidate columns for the stable admin-unit id, in preference order across sources/levels.
_ID_CANDIDATES = ("adm2_pcode", "adm1_pcode", "adm3_pcode", "shapeID", "GID_2", "GID_1", "pcode")
_ADM1_CANDIDATES = ("adm1_name", "NAME_1")


def fix_mojibake(text: object) -> object:
    """Repair the common latin-1/UTF-8 double-encoding seen in some shapefile names.

    e.g. ``"MboudÃ©"`` -> ``"Mboudé"``. Leaves already-correct strings unchanged.

    Args:
        text: A value that may be a mojibake string.

    Returns:
        The repaired string, or the input unchanged if it is not a fixable string.
    """
    if not isinstance(text, str) or "Ã" not in text:
        return text
    try:
        return text.encode("latin-1").decode("utf-8")
    except (UnicodeDecodeError, UnicodeEncodeError):
        return text


def normalize_country_nodes(iso3: str, gpkg_path: Path) -> gpd.GeoDataFrame:
    """Load one country GeoPackage and normalize it to the unified node schema.

    Args:
        iso3: Country code (injected as the ``iso3`` column; the source gpkgs do not all
            carry it).
        gpkg_path: Path to the country's ``*_admin*.gpkg``.

    Returns:
        GeoDataFrame with columns ``[iso3, adm2_name, adm2_id, adm1_name, population, lon,
        lat, geometry]`` in EPSG:4326.
    """
    gdf = gpd.read_file(gpkg_path)
    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(4326)

    adm2_id_col = next((c for c in _ID_CANDIDATES if c in gdf.columns), None)
    adm2_id = (
        gdf[adm2_id_col].astype(str)
        if adm2_id_col
        else pd.Series([f"{iso3}-{i}" for i in range(len(gdf))], index=gdf.index)
    )
    adm1_col = next((c for c in _ADM1_CANDIDATES if c in gdf.columns), None)
    adm1_name = gdf[adm1_col] if adm1_col else pd.Series([None] * len(gdf), index=gdf.index)

    name = gdf["name"] if "name" in gdf.columns else adm2_id
    pop = gdf["population"].fillna(0) if "population" in gdf.columns else pd.Series(0, index=gdf.index)

    # Geometric centroids for gravity distances. Centroid in a geographic CRS is slightly
    # off but fine at admin-2 scale; suppress the expected warning.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cent = gdf.geometry.centroid

    out = gpd.GeoDataFrame(
        {
            "iso3": iso3,
            "adm2_name": name.map(fix_mojibake).astype(str),
            "adm2_id": adm2_id.values,
            "adm1_name": adm1_name.map(fix_mojibake).values,
            "population": pop.astype(float).values,
            "lon": cent.x.values,
            "lat": cent.y.values,
            "geometry": gdf.geometry.values,
        },
        crs="EPSG:4326",
    )
    return out


def build_global_nodes(cfg: Config, isos: list[str] | None = None) -> gpd.GeoDataFrame:
    """Stack all available per-country node tables into one global GeoDataFrame.

    Args:
        cfg: Project configuration.
        isos: Optional explicit list of ISO3 codes; defaults to every country directory
            under ``data/countries/`` that has a GeoPackage.

    Returns:
        Global GeoDataFrame with a unique integer ``global_nodeid`` (0..N-1) plus the
        unified node schema, sorted by ``iso3`` then ``adm2_name``.

    Raises:
        FileNotFoundError: If no country GeoPackages are found at all.
    """
    if isos is None:
        isos = sorted(p.name for p in cfg.countries_dir.iterdir() if p.is_dir())

    frames: list[gpd.GeoDataFrame] = []
    for iso3 in isos:
        # Prefer admin-2; accept admin-1 fallback.
        gpkg = None
        for level in (2, 1):
            cand = cfg.countries_dir / iso3 / f"{iso3}_admin{level}.gpkg"
            if cand.exists():
                gpkg = cand
                break
        if gpkg is None:
            logger.info("nodes: no gpkg for %s, skipping", iso3)
            continue
        frames.append(normalize_country_nodes(iso3, gpkg))

    if not frames:
        raise FileNotFoundError(f"No country GeoPackages found under {cfg.countries_dir}")

    gdf = pd.concat(frames, ignore_index=True)
    gdf = gdf.sort_values(["iso3", "adm2_name"]).reset_index(drop=True)
    gdf.insert(0, "global_nodeid", range(len(gdf)))
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs="EPSG:4326")
    logger.info(
        "nodes: built %d global admin-2 nodes across %d countries (total pop=%s)",
        len(gdf), gdf["iso3"].nunique(), f"{gdf['population'].sum():,.0f}",
    )
    return gdf


def save_global_nodes(gdf: gpd.GeoDataFrame, cfg: Config) -> tuple[Path, Path]:
    """Write the global node table as GeoPackage (with geometry) and Parquet (attrs + WKB).

    Args:
        gdf: The global node GeoDataFrame.
        cfg: Project configuration.

    Returns:
        Tuple of (gpkg_path, parquet_path).
    """
    cfg.nodes_dir.mkdir(parents=True, exist_ok=True)
    gpkg_path = cfg.nodes_dir / "global_admin2_nodes.gpkg"
    parquet_path = cfg.nodes_dir / "global_admin2_nodes.parquet"
    gdf.to_file(gpkg_path, driver="GPKG")
    gdf.to_parquet(parquet_path)  # geopandas writes geometry as WKB in parquet
    logger.info("nodes: wrote %s and %s", gpkg_path.name, parquet_path.name)
    return gpkg_path, parquet_path


def load_global_nodes(cfg: Config) -> gpd.GeoDataFrame:
    """Load the previously-saved global node table from Parquet.

    Args:
        cfg: Project configuration.

    Returns:
        The global node GeoDataFrame.

    Raises:
        FileNotFoundError: If the node Parquet has not been built yet.
    """
    parquet_path = cfg.nodes_dir / "global_admin2_nodes.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(f"{parquet_path} not found; run build_global_nodes first")
    return gpd.read_parquet(parquet_path)
