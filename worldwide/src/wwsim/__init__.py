"""wwsim -- worldwide spatial data pipeline for LASER / Razer.

Turns open data into the inputs for a planet-scale metapopulation SEIR model:

1. **Acquire** admin-2 boundaries + matching 1km WorldPop population per UN member state,
   via the ``laser-init`` tool (:mod:`wwsim.acquire`).
2. **Nodes**: fuse the per-country GeoPackages into one global admin-2 node table
   (:mod:`wwsim.nodes`).
3. **Flights**: ingest open airport/route data (OurAirports + OpenFlights) and assign each
   airport to an admin-2 node (:mod:`wwsim.flights`, :mod:`wwsim.airports_assign`).
4. **Networks**: per-country gravity migration matrices (:mod:`wwsim.gravity`) and one
   global cross-border air-travel network over admin-2 nodes (:mod:`wwsim.air_network`),
   combined through a pluggable multi-modal layer (:mod:`wwsim.network`).
"""

from __future__ import annotations

from .config import Config, load_config
from .countries import UN_MEMBERS
from .logging import logger

__version__ = "0.1.0"

__all__ = ["Config", "load_config", "UN_MEMBERS", "logger", "__version__"]
