#!/usr/bin/env python3
"""Convert the EnglandWalesMeasles.py dataset into a shareable CSV.

The source module (`EnglandWalesMeasles.py`) stores 954 England & Wales
registration districts ("places") as Python `Place` objects, each holding
per-year time series (population, births, cases over 1944-1964) plus a fixed
latitude/longitude. That Python-only representation is awkward to consume from
R, so this script flattens the static per-place attributes into a CSV with one
row per node:

    name, population, latitude, longitude

`population` is the first-year (1944) population, used as each node's initial
population for the SIR example. The full time series remain available in the
source module if richer scenarios are needed later.

Run from the package root:  python3 examples/data/convert_measles.py
"""

import csv
import importlib.util
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "EnglandWalesMeasles.py"
DEST = HERE / "EnglandWalesMeasles_places.csv"


def load_data():
    """Import the EnglandWalesMeasles module from its file path.

    Returns:
        The module's `data` container (with `.placenames` and `.places`).
    """
    logger.info("Loading dataset from %s", SOURCE)
    spec = importlib.util.spec_from_file_location("englandwalesmeasles", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.data


def main():
    data = load_data()
    logger.info("Found %d places", len(data.placenames))

    with DEST.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "population", "latitude", "longitude"])
        for name in data.placenames:
            place = data.places[name]
            # population is a per-year array; take year 0 (1944) as the initial value.
            writer.writerow([
                name,
                int(place.population[0]),
                float(place.latitude),
                float(place.longitude),
            ])

    logger.info("Wrote %s", DEST)


if __name__ == "__main__":
    main()
