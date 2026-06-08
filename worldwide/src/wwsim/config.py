"""Central configuration: filesystem layout and default model parameters.

All paths are anchored at the project root (the ``worldwide/`` directory, the parent of
``src/``) so scripts run from anywhere produce a consistent on-disk layout:

```
worldwide/
  data/
    countries/<ISO3>/<ISO3>_admin2.gpkg   # per-country admin-2 + population (from laser-init)
    flights/                              # OurAirports / OpenFlights raw downloads
    interim/                              # cached intermediate artifacts
  output/
    nodes/        global_admin2_nodes.*   # unified worldwide admin-2 node table
    networks/     <ISO3>_gravity.*        # per-country intra-country migration matrices
                  global_air_network.*    # global cross-border air network (admin-2 nodes)
                  global_combined_network.* # intra (block-diagonal) + inter (air) combined
    plots/                                # validation figures
```

Defaults live here as a single :class:`Config` dataclass so scripts and tests share one
source of truth. Anything here can be overridden from a YAML file via :func:`load_config`.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .logging import logger

# Project root = parent of the directory holding this file's package (src/wwsim -> src -> root).
PROJECT_ROOT = Path(__file__).resolve().parents[2]


@dataclass
class GravityParams:
    """Parameters for the gravity migration model.

    The connection weight from source node ``i`` to destination node ``j`` is

        ``w_ij = k * P_i^a * P_j^b / D_ij^c``

    where ``P`` is population and ``D`` is great-circle distance (km). Defaults match
    laser-init's generated models (``k=500, a=1, b=1, c=2``).

    Attributes:
        k: Overall scaling constant.
        a: Exponent on the source population.
        b: Exponent on the destination population.
        c: Exponent on distance (decay).
        min_distance_km: Distance floor (km) to avoid division blow-up for near-coincident
            centroids and self-distance.
    """

    k: float = 500.0
    a: float = 1.0
    b: float = 1.0
    c: float = 2.0
    min_distance_km: float = 1.0


@dataclass
class AirNetworkParams:
    """Parameters for the global cross-border air-travel network.

    Attributes:
        top_n_airports: Keep only the busiest N airports by passenger-volume proxy before
            building edges. ``None`` keeps all airports.
        passenger_proxy: How to weight an airport->airport route. ``"seats"`` sums modelled
            seats per carrier-route (uses the aircraft seat table); ``"routes"`` counts
            distinct carrier-routes (route multiplicity).
        exclude_codeshares: Drop OpenFlights rows flagged as codeshares to avoid counting
            one physical flight multiple times.
        direct_only: Keep only non-stop routes (``stops == 0``).
        intra_country_air: If ``False`` (the goal), drop edges whose endpoints are in the
            same country -- the global network carries only cross-border air travel.
    """

    top_n_airports: int | None = 1000
    passenger_proxy: str = "seats"
    exclude_codeshares: bool = True
    direct_only: bool = True
    intra_country_air: bool = False


@dataclass
class Config:
    """Top-level configuration bundle.

    Attributes:
        project_root: Project root directory.
        year: Target data year (population raster + demographics).
        adm_level: Preferred administrative level (2); the acquirer falls back to 1.
        shape_source_order: Waterfall of shapefile sources tried per country.
        gravity: Gravity-model parameters.
        air: Air-network parameters.
    """

    project_root: Path = PROJECT_ROOT
    year: int = 2015
    adm_level: int = 2
    shape_source_order: tuple[str, ...] = ("unocha", "geoboundaries", "gadm")
    gravity: GravityParams = field(default_factory=GravityParams)
    air: AirNetworkParams = field(default_factory=AirNetworkParams)

    # ---- derived paths (properties so they always reflect project_root) ----
    @property
    def data_dir(self) -> Path:
        return self.project_root / "data"

    @property
    def countries_dir(self) -> Path:
        """Per-country laser-init output (one subdir per ISO3)."""
        return self.data_dir / "countries"

    @property
    def flights_dir(self) -> Path:
        return self.data_dir / "flights"

    @property
    def interim_dir(self) -> Path:
        return self.data_dir / "interim"

    @property
    def output_dir(self) -> Path:
        return self.project_root / "output"

    @property
    def nodes_dir(self) -> Path:
        return self.output_dir / "nodes"

    @property
    def networks_dir(self) -> Path:
        return self.output_dir / "networks"

    @property
    def plots_dir(self) -> Path:
        return self.output_dir / "plots"

    def country_gpkg(self, iso3: str) -> Path:
        """Path to a country's admin-2 GeoPackage (may be admin1 if admin2 unavailable).

        Args:
            iso3: ISO 3166-1 alpha-3 country code.

        Returns:
            Expected GeoPackage path under ``data/countries/<ISO3>/``.
        """
        return self.countries_dir / iso3 / f"{iso3}_admin{self.adm_level}.gpkg"

    def ensure_dirs(self) -> None:
        """Create all output/data directories if missing (idempotent)."""
        for d in (
            self.data_dir,
            self.countries_dir,
            self.flights_dir,
            self.interim_dir,
            self.output_dir,
            self.nodes_dir,
            self.networks_dir,
            self.plots_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)


def load_config(path: Path | str | None = None) -> Config:
    """Load configuration, overlaying values from a YAML file when provided.

    Args:
        path: Optional path to a YAML file with any subset of top-level keys
            (``year``, ``adm_level``, ``shape_source_order``, ``gravity``, ``air``).

    Returns:
        A :class:`Config`. When ``path`` is ``None`` (or missing) the built-in defaults
        are returned unchanged.

    Raises:
        FileNotFoundError: If ``path`` is given but does not exist.
    """
    cfg = Config()
    if path is None:
        return cfg

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    raw = yaml.safe_load(path.read_text()) or {}
    for key in ("year", "adm_level"):
        if key in raw:
            setattr(cfg, key, raw[key])
    if "shape_source_order" in raw:
        cfg.shape_source_order = tuple(raw["shape_source_order"])
    if "gravity" in raw:
        cfg.gravity = GravityParams(**{**dataclasses.asdict(cfg.gravity), **raw["gravity"]})
    if "air" in raw:
        cfg.air = AirNetworkParams(**{**dataclasses.asdict(cfg.air), **raw["air"]})

    logger.info("Loaded config overrides from %s", path)
    return cfg
