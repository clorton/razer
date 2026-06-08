# worldwide — planet-scale spatial inputs for LASER / Razer

`wwsim` builds the spatial inputs for a **worldwide SEIR metapopulation simulation** of
COVID-19 (or any directly-transmitted pathogen) from **open data, vintage ~2015**:

- **Admin-2 nodes** for all **193 UN member states** — boundaries + matching **1 km
  WorldPop 2015** population, aggregated per admin-2 unit.
- **Per-country intra-country migration matrices** (one gravity network per country).
- **One global cross-border air-travel network**: airports link **admin-2 nodes across
  country borders**, limited to the **top-N airports by passenger volume**, with all
  airports in the same admin-2 unit aggregated into one node.
- A **pluggable multi-modal layer** so passenger **rail** (and a licensed **OAG** O-D
  feed) can be summed into the inter-nation network later without touching downstream code.

## Data sources (all open)

| Layer | Source | Via | License |
|---|---|---|---|
| Admin-2 boundaries | UNOCHA COD-AB → geoBoundaries → GADM (waterfall) | [`laser-init`](https://github.com/laser-base/laser-init) | per-source (COD per-country; geoBoundaries CC-BY) |
| Population raster | WorldPop Global 1 km, 2015, UN-adjusted constrained | `laser-init` | CC-BY 4.0 |
| Zonal aggregation | [RasterToolkit](https://github.com/InstituteforDiseaseModeling/RasterToolkit) | `laser-init` | MIT |
| Airports (coords, ISO country) | [OurAirports](https://ourairports.com/data/) | `wwsim.flights` | Public domain |
| Air routes (O-D, ~2014) | [OpenFlights](https://github.com/jpatokal/openflights) | `wwsim.flights` | ODbL |
| Demographics (optional) | UN World Population Prospects 2024 | `laser-init` | CC-BY 3.0 IGO |

> **Why not OAG?** OAG passenger origin-destination volumes are a commercial product; no
> free equivalent exists. We substitute an OpenFlights route network weighted by a
> **seat-capacity passenger proxy**, behind an adapter so a licensed OAG file drops in
> later (`wwsim.oag`).

## Install

Requires [`uv`](https://docs.astral.sh/uv/).

```shell
# from worldwide/
uv venv .venv --python 3.12
uv pip install git+https://github.com/laser-base/laser-init   # the acquisition tool (+ rastertoolkit)
uv pip install -e .                      # this package (wwsim)
```

`laser-init` is cloned into `vendor/` (see repo setup); it is not on PyPI.

## Pipeline

Each step is a script under `scripts/`, and each wraps a `wwsim` module:

```shell
.venv/bin/python scripts/01_acquire_countries.py        # admin-2 gpkgs for 193 countries
.venv/bin/python scripts/02_build_nodes.py              # unified global admin-2 node table
.venv/bin/python scripts/03_ingest_flights.py           # OurAirports + OpenFlights
.venv/bin/python scripts/04_assign_airports.py          # airport -> admin-2 (point-in-polygon)
.venv/bin/python scripts/05_intra_country_matrices.py   # per-country gravity networks
.venv/bin/python scripts/06_global_air_network.py       # global cross-border air network
.venv/bin/python scripts/07_combine_network.py          # intra + inter combined global matrix
```

See **`USAGE.md`** for the full step-by-step guide (flags, top-N selection, configuration,
loading into LASER/Razer), `docs/datasources.md` for source details, and `docs/pipeline.md`
for the data model.

## Outputs

- `output/nodes/global_admin2_nodes.{parquet,gpkg}` — every admin-2 node worldwide:
  `global_nodeid, iso3, adm2_name, adm2_id, population, lon, lat, geometry`.
- `output/networks/<ISO3>_gravity.{parquet,npz}` — per-country intra-country matrices.
- `output/networks/global_air_network.{parquet,npz}` — global cross-border air network.
- `output/networks/global_combined_network.npz` — intra (block-diagonal) + inter (air).
- `output/plots/` — choropleths, network maps, degree distributions.
