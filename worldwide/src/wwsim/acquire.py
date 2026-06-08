"""Acquire admin-2 boundaries + matching 1km WorldPop population per country.

This module drives `laser-init` (the laser-base bootstrap tool) over every UN member
state. For each country it produces a GeoPackage ``<ISO3>_admin2.gpkg`` whose rows are
admin-2 units carrying a ``population`` column aggregated from the 2015 1km WorldPop
raster by RasterToolkit.

**Source waterfall.** Per the configured order (default ``unocha -> geoboundaries ->
gadm``) and level fallback (admin-2, then admin-1), the first source that yields features
wins. This guarantees coverage for all 193 members: UNOCHA COD-AB is authoritative but
only reaches admin-2 for roughly half the world, geoBoundaries fills the rest, and the
admin-1 fallback covers the ~dozen microstates with no admin-2 anywhere.

**UNOCHA efficiency.** laser-init's stock UNOCHA transformer re-reads the 1-2 GB global
geodatabase on every call. Tried first for 193 countries that would mean ~193 full reads.
Instead :class:`Acquirer` downloads and reads each global admin-level layer **once**, caches
the GeoDataFrame in memory, and filters per country -- turning the per-country UNOCHA step
into an in-memory selection plus a small raster clip.

The WorldPop raster is downloaded **once per country** (it is source-independent) and
reused across whichever shape source ends up winning.
"""

from __future__ import annotations

import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import geopandas as gpd

# Reuse laser-init's downloaders/transformers so we genuinely "invoke" the tool.
from laser.init.cli import download_raster_data, download_shape_data, transform_shape_and_raster_data
from laser.init.extractors.unocha import UnochaExtractor
from laser.init.transformers.unocha import read_gdb_quietly
from laser.init.utils import clip_quietly

from .config import Config
from .logging import logger
from .zonal import zonal_sum


class CountryNotAvailable(Exception):
    """Raised when a (source, level) combination has no features for a country."""


@dataclass
class AcquireResult:
    """Outcome of acquiring one country.

    Attributes:
        iso3: Country code.
        status: ``"ok"``, ``"skipped"`` (already present), or ``"failed"``.
        source: Winning shape source (``unocha``/``geoboundaries``/``gadm``), or ``None``.
        level: Administrative level actually used (2 or 1), or ``None``.
        gpkg: Path to the produced GeoPackage, or ``None``.
        n_units: Number of admin units in the GeoPackage, or ``None``.
        population: Total summed population, or ``None``.
        error: Error summary if failed.
    """

    iso3: str
    status: str
    source: str | None = None
    level: int | None = None
    gpkg: Path | None = None
    n_units: int | None = None
    population: float | None = None
    error: str | None = None


@dataclass
class Acquirer:
    """Stateful acquirer that caches the UNOCHA global layers across countries.

    Args:
        cfg: Project configuration (paths, year, source order).
        levels: Administrative levels to try, in order (default admin-2 then admin-1).
    """

    cfg: Config
    levels: tuple[int, ...] = (2, 1)
    _unocha_layers: dict[int, gpd.GeoDataFrame] = field(default_factory=dict, repr=False)
    _unocha_gdb_dir: Path | None = field(default=None, repr=False)

    # ------------------------------------------------------------------ UNOCHA
    def _ensure_unocha_gdb(self) -> Path:
        """Download (once) and extract (once) the UNOCHA global geodatabase.

        Returns:
            Path to the extracted ``.gdb`` directory.
        """
        if self._unocha_gdb_dir is not None:
            return self._unocha_gdb_dir

        zip_path = UnochaExtractor().extract(self.cfg.year, self.cfg.adm_level, self.cfg.year)
        gdb_dir = zip_path.parent / zip_path.stem
        if not gdb_dir.exists():
            logger.info("Extracting UNOCHA global geodatabase %s (one time)...", zip_path.name)
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(path=zip_path.parent)
        self._unocha_gdb_dir = gdb_dir
        return gdb_dir

    def _unocha_layer(self, level: int) -> gpd.GeoDataFrame:
        """Return the global admin-``level`` layer, reading it from disk at most once.

        Args:
            level: Administrative level (0-4).

        Returns:
            The full global GeoDataFrame for that level.
        """
        if level not in self._unocha_layers:
            gdb_dir = self._ensure_unocha_gdb()
            logger.info("Reading UNOCHA global admin%d layer (one time)...", level)
            self._unocha_layers[level] = read_gdb_quietly(gdb_dir, layer_name=f"admin{level}")
        return self._unocha_layers[level]

    def _build_unocha(self, iso3: str, level: int, raster: Path, out_dir: Path) -> Path:
        """Build a country GeoPackage from the preloaded UNOCHA global layer.

        Mirrors laser-init's ``UnochaTransformer`` but reuses the cached global frame.

        Args:
            iso3: Country code.
            level: Administrative level.
            raster: Path to the country's WorldPop raster.
            out_dir: Output directory for the GeoPackage.

        Returns:
            Path to the written GeoPackage.

        Raises:
            CountryNotAvailable: If the global layer has no rows for ``iso3``.
        """
        layer = self._unocha_layer(level)
        country = layer[layer.iso3 == iso3]
        if len(country) == 0:
            raise CountryNotAvailable(f"UNOCHA admin{level} has no features for {iso3}")

        names = [f"adm{i}_name" for i in range(level + 1)]
        pcode = f"adm{level}_pcode"
        country = country[names + [pcode, "geometry"]].copy()
        country["nodeid"] = list(range(len(country)))
        country["name"] = country[f"adm{level}_name"]

        # RasterToolkit needs a shapefile keyed by a stable id (the p-code).
        with tempfile.TemporaryDirectory() as tmp:
            shp = Path(tmp) / f"{iso3}_admin{level}.shp"
            country.to_file(shp, driver="ESRI Shapefile", engine="pyogrio")
            pop = clip_quietly(raster, shp, shape_attr=pcode)
        country["population"] = country[pcode].map(pop)

        out = out_dir / f"{iso3}_admin{level}.gpkg"
        country.to_file(out, driver="GPKG")
        logger.info("UNOCHA: wrote %s (%d admin%d units)", out.name, len(country), level)
        return out

    # ------------------------------------------------------- geoBoundaries/GADM
    def _build_via_laserinit(
        self, source: str, iso3: str, level: int, raster: Path, out_dir: Path
    ) -> Path:
        """Build a country GeoPackage via laser-init's per-country extractor + transformer.

        Args:
            source: ``"geoboundaries"`` or ``"gadm"``.
            iso3: Country code.
            level: Administrative level.
            raster: Path to the country's WorldPop raster.
            out_dir: Output directory.

        Returns:
            Path to the written GeoPackage.
        """
        shape = download_shape_data(iso3, level, self.cfg.year, source)
        return transform_shape_and_raster_data(source, shape, iso3, level, raster, out_dir)

    # ----------------------------------------- rasterio fallback (antimeridian etc.)
    def _shapes_for(self, source: str, iso3: str, level: int) -> tuple[gpd.GeoDataFrame, str]:
        """Load a country's boundary polygons (no clipping) for a given source/level.

        Args:
            source: ``"unocha"``, ``"geoboundaries"`` or ``"gadm"``.
            iso3: Country code.
            level: Administrative level.

        Returns:
            Tuple of (GeoDataFrame in EPSG:4326 with at least ``name`` + an id column +
            ``geometry``, name of the id column).

        Raises:
            CountryNotAvailable: If UNOCHA has no features for the country/level.
        """
        if source == "unocha":
            layer = self._unocha_layer(level)
            c = layer[layer.iso3 == iso3]
            if len(c) == 0:
                raise CountryNotAvailable(f"UNOCHA admin{level} has no features for {iso3}")
            names = [f"adm{i}_name" for i in range(level + 1)]
            pcode = f"adm{level}_pcode"
            g = gpd.GeoDataFrame(c[names + [pcode, "geometry"]].copy(), geometry="geometry", crs=c.crs)
            g["name"] = g[f"adm{level}_name"]
            return _to_4326(g), pcode

        if source == "geoboundaries":
            from laser.init.extractors.geoboundaries import GeoBoundariesExtractor

            path = GeoBoundariesExtractor().extract(iso3, level, self.cfg.year)
            g = gpd.read_file(path, layer=f"geoBoundaries-{iso3.upper()}-ADM{level}")
            g = g[["shapeName", "shapeID", "geometry"]].copy()
            g["name"] = g["shapeName"]
            return _to_4326(g), "shapeID"

        if source == "gadm":
            from laser.init.extractors.gadm import GadmExtractor

            path = GadmExtractor().extract(iso3, level, self.cfg.year)
            g = gpd.read_file(path, layer=f"gadm41_{iso3.upper()}_{level}")
            gid, nm = f"GID_{level}", f"NAME_{level}"
            cols = [c for c in (nm, gid) if c in g.columns] + ["geometry"]
            g = g[cols].copy()
            g["name"] = g[nm] if nm in g.columns else g[gid]
            return _to_4326(g), gid

        raise ValueError(f"unknown source {source!r}")

    def _build_rasterio(self, iso3: str, level: int, source: str, raster: Path, out_dir: Path) -> Path:
        """Build a country GeoPackage using rasterio zonal sums (no RasterToolkit).

        Used for countries whose WorldPop file is delivered on a global canvas (tie point
        ``x0 = -180``), which RasterToolkit's strict assertion rejects.

        Args:
            iso3: Country code.
            level: Administrative level.
            source: Shape source.
            raster: Path to the WorldPop raster.
            out_dir: Output directory.

        Returns:
            Path to the written GeoPackage.
        """
        gdf, id_col = self._shapes_for(source, iso3, level)
        pop = zonal_sum(raster, gdf, id_col)
        gdf["population"] = gdf[id_col].map(pop)
        out = out_dir / f"{iso3}_admin{level}.gpkg"
        gdf.to_file(out, driver="GPKG")
        logger.info("rasterio: wrote %s (%d units via %s)", out.name, len(gdf), source)
        return out

    def acquire_country_rasterio(self, iso3: str, force: bool = True) -> AcquireResult:
        """Acquire one country using the rasterio zonal-sum fallback for population.

        Same source/level waterfall as :meth:`acquire_country`, but population is summed
        with rasterio so antimeridian / global-canvas rasters work.

        Args:
            iso3: Country code.
            force: Rebuild even if a gpkg exists (default True -- this is the repair path).

        Returns:
            An :class:`AcquireResult`.
        """
        out_dir = self.cfg.countries_dir / iso3
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            raster = self._download_raster(iso3)
        except Exception as exc:  # noqa: BLE001
            return AcquireResult(iso3, "failed", error=f"raster: {exc}")

        last_error: str | None = None
        for source in self.cfg.shape_source_order:
            for level in self.levels:
                try:
                    gpkg = self._build_rasterio(iso3, level, source, raster, out_dir)
                    n, pop = _summarize_gpkg(gpkg)
                    # Some R2025A country files are empty global canvases (e.g. KIR).
                    # If the primary raster yields no population, fall back to the older
                    # 2000-2020 product and rebuild this same source/level.
                    if pop == 0:
                        logger.info("%s: primary raster empty -> 2000-2020 fallback", iso3)
                        alt = self._download_raster(iso3, "global2020")
                        gpkg = self._build_rasterio(iso3, level, source, alt, out_dir)
                        n, pop = _summarize_gpkg(gpkg)
                    logger.info("%s: OK(rasterio) via %s admin%d (%d units, pop=%s)",
                                iso3, source, level, n, f"{pop:,.0f}")
                    return AcquireResult(iso3, "ok", source, level, gpkg, n, pop)
                except Exception as exc:  # noqa: BLE001
                    last_error = f"{source}/admin{level}: {exc}"
                    logger.info("%s: rasterio %s -- trying next", iso3, last_error)
        return AcquireResult(iso3, "failed", error=last_error)

    # ----------------------------------------------------------------- raster
    def _download_raster(self, iso3: str, product: str = "r2025a") -> Path:
        """Download the country's 2015 1km WorldPop raster (cached by laser-init).

        Args:
            iso3: Country code.
            product: ``"r2025a"`` for the Global 2015-2030 R2025A constrained UN-adjusted
                product (laser-init's default), or ``"global2020"`` for the older Global
                2000-2020 1km UN-adjusted product -- a fallback for the handful of countries
                whose R2025A file is an empty global canvas (e.g. Kiribati).

        Returns:
            Path to the GeoTIFF.

        Raises:
            RuntimeError: If the raster cannot be downloaded.
        """
        if product == "r2025a":
            return download_raster_data(iso3, self.cfg.year, "worldpop")
        if product == "global2020":
            from laser.init.config import configuration as licfg
            from laser.init.config import default_cache_directory
            from laser.init.utils import download_file

            cache_root = Path(licfg.get("cache_dir", default_cache_directory))
            fn = f"{iso3.lower()}_ppp_{self.cfg.year}_1km_Aggregated_UNadj.tif"
            url = (
                "https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/"
                f"{self.cfg.year}/{iso3}/{fn}"
            )
            return download_file(url, cache_dir=cache_root, dest_dir=Path("WorldPop"))
        raise ValueError(f"unknown raster product {product!r}")

    # --------------------------------------------------------------- per country
    def acquire_country(self, iso3: str, force: bool = False) -> AcquireResult:
        """Acquire one country with the full source/level waterfall.

        Args:
            iso3: ISO 3166-1 alpha-3 code.
            force: Re-acquire even if a GeoPackage already exists.

        Returns:
            An :class:`AcquireResult` describing the outcome.
        """
        out_dir = self.cfg.countries_dir / iso3
        out_dir.mkdir(parents=True, exist_ok=True)

        # Resumability: accept an existing gpkg at any tried level.
        if not force:
            for level in self.levels:
                existing = out_dir / f"{iso3}_admin{level}.gpkg"
                if existing.exists():
                    n, pop = _summarize_gpkg(existing)
                    logger.info("%s: already present (%s), skipping", iso3, existing.name)
                    return AcquireResult(iso3, "skipped", None, level, existing, n, pop)

        try:
            raster = self._download_raster(iso3)
        except Exception as exc:  # noqa: BLE001 - record and bail; raster is prerequisite
            logger.warning("%s: raster download failed: %s", iso3, exc)
            return AcquireResult(iso3, "failed", error=f"raster: {exc}")

        last_error: str | None = None
        for source in self.cfg.shape_source_order:
            for level in self.levels:
                try:
                    if source == "unocha":
                        gpkg = self._build_unocha(iso3, level, raster, out_dir)
                    else:
                        gpkg = self._build_via_laserinit(source, iso3, level, raster, out_dir)
                    n, pop = _summarize_gpkg(gpkg)
                    logger.info(
                        "%s: OK via %s admin%d (%d units, pop=%s)", iso3, source, level, n, f"{pop:,.0f}"
                    )
                    return AcquireResult(iso3, "ok", source, level, gpkg, n, pop)
                except Exception as exc:  # noqa: BLE001 - try the next source/level
                    last_error = f"{source}/admin{level}: {exc}"
                    logger.info("%s: %s -- trying next", iso3, last_error)

        logger.warning("%s: FAILED all sources/levels", iso3)
        return AcquireResult(iso3, "failed", error=last_error)


def _to_4326(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Return the GeoDataFrame in EPSG:4326 (reproject only if needed)."""
    if gdf.crs is None:
        return gdf.set_crs("EPSG:4326")
    if gdf.crs.to_epsg() != 4326:
        return gdf.to_crs(4326)
    return gdf


def _summarize_gpkg(path: Path) -> tuple[int, float]:
    """Return (row count, total population) for a country GeoPackage.

    Args:
        path: Path to a ``*_admin*.gpkg``.

    Returns:
        Tuple of (number of admin units, summed population). Population is ``0.0`` if the
        column is missing or all-null.
    """
    gdf = gpd.read_file(path)
    pop = float(gdf["population"].fillna(0).sum()) if "population" in gdf.columns else 0.0
    return len(gdf), pop
