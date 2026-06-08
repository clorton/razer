"""Rasterio-based zonal population sums -- a robust fallback to RasterToolkit.

RasterToolkit is fast and dependency-light, but it asserts the raster tie point lies
strictly inside ``-180 < x0 < 180`` / ``-85 < y0 < 85``. A handful of WorldPop R2025A
country files are delivered on a **full global canvas** (43200 px wide, tie point exactly
``x0 = -180.0``) -- typically antimeridian-crossing countries (Russia, Fiji, Kiribati,
Tuvalu) and a few packaging quirks (Eritrea). Those fail the strict assertion even though
the data is correct.

This module sums population per polygon directly with rasterio, reading **only each
polygon's window** (so even a global-canvas raster is cheap), and is used as a fallback by
the rebuild path for those countries.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.mask import mask as rio_mask

from .logging import logger


def zonal_sum(raster_path: Path, gdf: gpd.GeoDataFrame, id_col: str) -> dict[str, float]:
    """Sum raster values within each polygon, keyed by an id column.

    Args:
        raster_path: Path to the population GeoTIFF (assumed EPSG:4326, persons/pixel).
        gdf: GeoDataFrame of polygons (EPSG:4326) with a unique ``id_col``.
        id_col: Column whose values key the returned dict.

    Returns:
        Mapping ``id -> summed population`` (rounded to int-like float). Polygons that fall
        entirely outside the raster get ``0.0``. Nodata and implausible negative sentinels
        (e.g. ``-99999``, ``-3.4e38``) are excluded; only non-negative cells are summed.
    """
    out: dict[str, float] = {}
    with rasterio.open(raster_path) as src:
        nodata = src.nodata
        for _, row in gdf.iterrows():
            geom = [row.geometry.__geo_interface__]
            try:
                arr, _ = rio_mask(src, geom, crop=True, all_touched=False, filled=True)
            except ValueError:
                # "Input shapes do not overlap raster" -> polygon outside coverage.
                out[row[id_col]] = 0.0
                continue
            a = arr[0].astype("float64").ravel()
            if nodata is not None:
                a = a[a != nodata]
            # Population is non-negative; this also drops huge-negative nodata sentinels.
            a = a[(a >= 0) & np.isfinite(a)]
            out[row[id_col]] = float(np.round(a.sum()))
    logger.info("zonal_sum: summed %d polygons from %s", len(out), Path(raster_path).name)
    return out
