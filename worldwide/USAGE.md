# USAGE — `wwsim` worldwide pipeline

A practical, end-to-end guide. For *what* the data is and its licenses see `README.md` and
`docs/datasources.md`; for the data model see `docs/pipeline.md`.

---

## 1. Install

Requires [`uv`](https://docs.astral.sh/uv/). From `worldwide/`:

```shell
uv venv .venv --python 3.12
uv pip install git+https://github.com/laser-base/laser-init    # acquisition tool (+ rastertoolkit, geopandas, …)
uv pip install -e .                       # this package (wwsim)
uv pip install rasterio                   # zonal-sum fallback for global-canvas rasters
```

`laser-init` is cloned under `vendor/` (it is not on PyPI). Downloaded source data is cached
in `~/.laser/cache/` (outside the repo) and reused across runs.

---

## 2. Run everything

```shell
# Full build from scratch (acquire 193 countries, ingest flights, build all networks, plot):
.venv/bin/python scripts/run_all.py --with-acquire --with-flights

# If countries + flights are already acquired, just (re)build nodes → networks → plots:
.venv/bin/python scripts/run_all.py
```

`run_all.py` runs steps 02→08 (and 01+01b, 03 with the flags). Every step is **resumable**:
acquisition skips countries whose GeoPackage already exists; downloads are cached.

---

## 3. Run step by step

| # | Script | Does | Key outputs |
|---|---|---|---|
| 01 | `01_acquire_countries.py` | Admin-2 boundaries + 1 km WorldPop 2015 per country (via laser-init) | `data/countries/<ISO3>/<ISO3>_admin2.gpkg`, `output/nodes/acquisition_manifest.csv` |
| 01b | `01b_fix_failed.py` | Repair failures with the rasterio zonal-sum fallback | rewrites the failed gpkgs |
| 02 | `02_build_nodes.py` | Fuse per-country gpkgs into the global node table | `output/nodes/global_admin2_nodes.{parquet,gpkg}` |
| 03 | `03_ingest_flights.py` | Download OurAirports + OpenFlights, build airport edges | `data/interim/airport_edges.parquet`, `airports_master.parquet` |
| 04 | `04_assign_airports.py` | Airport → admin-2 node (point-in-polygon) | `data/interim/airport_node_assignment.parquet` |
| 05 | `05_intra_country_matrices.py` | Per-country gravity matrices | `output/networks/<ISO3>_gravity.{npz,parquet}` |
| 06 | `06_global_air_network.py` | Global cross-border air network | `output/networks/global_air_network.{parquet,npz}` |
| 07 | `07_combine_network.py` | Intra + inter combined matrix | `output/networks/global_combined_network.npz` |
| 08 | `08_plots.py` | Validation figures | `output/plots/*.png` |

### Useful invocations

```shell
# Acquire only some countries, or force a fresh build, or override the source waterfall:
.venv/bin/python scripts/01_acquire_countries.py NGA KEN FRA
.venv/bin/python scripts/01_acquire_countries.py COM --force
.venv/bin/python scripts/01_acquire_countries.py --source-order geoboundaries gadm

# Repair specific countries with the rasterio fallback:
.venv/bin/python scripts/01b_fix_failed.py RUS FJI KIR

# Per-country gravity choropleth for a chosen country:
.venv/bin/python scripts/08_plots.py --country IND
```

---

## 4. Choosing the top-N busiest airports

The inter-country network is built from the **top-N airports by passenger-volume proxy**
(`Config.air.top_n_airports`, default 1000). This affects **only** the air network (and the
combined matrix); the node table and the 193 gravity matrices are untouched.

```shell
# One air network at a chosen cut-off (canonical filename):
.venv/bin/python scripts/06_global_air_network.py --top-n 250

# Several cut-offs at once → suffixed files global_air_network_top{N}.{parquet,npz}:
.venv/bin/python scripts/06_global_air_network.py --top-n 1000 500 250 100

# No cap (keep all airports):
.venv/bin/python scripts/06_global_air_network.py --top-n 0

# Combined matrices per cut-off (rebuilds the air layer in-process; step 06 not required):
.venv/bin/python scripts/07_combine_network.py --top-n 1000 500 250 100
#   → global_combined_network_top{N}.npz
# A single value (or none) writes the canonical global_combined_network.npz.
```

Observed concentration (first full build): the busiest **500** airports already carry ~88%
of cross-border seat volume; **250** ≈ 63%; **100** ≈ 30%.

---

## 5. Configuration

Defaults live in `src/wwsim/config.py` (`Config`, `GravityParams`, `AirNetworkParams`).
Override via a YAML file passed to scripts that accept `--config`:

```yaml
# my_config.yaml
year: 2015
adm_level: 2
shape_source_order: [unocha, geoboundaries, gadm]
gravity:
  k: 500      # w_ij = k * P_i^a * P_j^b / D_ij^c
  a: 1
  b: 1
  c: 2
air:
  top_n_airports: 500
  passenger_proxy: seats     # or "routes" (route multiplicity)
  exclude_codeshares: true
  direct_only: true
  intra_country_air: false   # keep only between-country air edges
```

```shell
.venv/bin/python scripts/01_acquire_countries.py --config my_config.yaml
```

---

## 6. Outputs reference

- `output/nodes/global_admin2_nodes.parquet` (and `.gpkg`) — the **node table**, one row per
  admin-2 unit: `global_nodeid, iso3, adm2_name, adm2_id, adm1_name, population, lon, lat,
  geometry`. `global_nodeid` (0…N-1) is the row/column index of every matrix.
- `output/networks/<ISO3>_gravity.npz` — per-country gravity: `ids` (global node ids) +
  `matrix` (dense float32). `<ISO3>_gravity.parquet` — the same as an edge list.
- `output/networks/global_air_network.{parquet,npz}` — cross-border admin-2 air edges
  (`src_global_nodeid, dst_global_nodeid, src_iso3, dst_iso3, weight, n_airport_pairs,
  n_carriers`) and the sparse matrix.
- `output/networks/global_combined_network.npz` — combined `(N, N)` sparse matrix
  (gravity block-diagonal + air off-diagonal).
- `output/nodes/acquisition_manifest.csv` — per country: status, source, level, units, pop.
- `output/plots/*.png` — choropleths and network maps.

---

## 7. Use it in LASER / Razer

```python
import geopandas as gpd, scipy.sparse as sp
nodes = gpd.read_parquet("output/nodes/global_admin2_nodes.parquet")  # per-patch pop + coords
M = sp.load_npz("output/networks/global_combined_network.npz")        # (N, N), N == len(nodes)
pop = nodes["population"].to_numpy()
```

`M[i, j]` is the directed coupling node *i* → node *j*. **Scaling caveat:** the gravity layer
(raw `k·Pᵃ·Pᵇ/Dᶜ`) and the air layer (seat proxy) are in different units. Out of the box
both `*_scale = 1.0`, so gravity dominates; for a single mixing kernel, set
`--gravity-scale`/`--air-scale` in step 07 (e.g. row-normalize each layer, or boost air), or
consume the two layers separately. `wwsim.gravity.row_normalize` converts a layer to
out-migration probabilities.

---

## 8. Swapping in real data

- **OAG (licensed) air O-D** — drop your CSV (`origin, destination, passengers`) through
  `wwsim.oag.load_oag_airport_edges(path, airports)`; it yields the same airport-edge schema
  `build_airport_edges` produces, so step 06 consumes it unchanged.
- **Rail** — `wwsim.rail.rail_mode_from_csv(path, nodes)` accepts either a node-id edge CSV
  or a station-coordinate CSV (snapped to admin-2), returning a `ModeNetwork` to add via
  `--rail-csv` in step 07. See the module docstring for candidate open sources.

---

## 9. Tests

```shell
.venv/bin/python -m pytest tests/ -q        # 33 given-when-then tests
```

---

## 10. Known caveats

- **Passenger proxy, not measured pax** — no free OAG equivalent; weights are relative
  seat-capacity. Calibrate against World Bank `IS.AIR.PSGR` for absolute volumes.
- **Constrained WorldPop** slightly undercounts (e.g. RUS 127.7M vs ~144M). Switch product
  in `wwsim.acquire.Acquirer._download_raster` if you prefer unconstrained.
- **Antimeridian / empty-canvas rasters** (RUS, FJI, KIR, TUV, ERI) are handled by the
  rasterio fallback (`01b_fix_failed.py`); KIR uses the 2000-2020 raster product.
- **220/9,056 airports unassigned** — mostly in dependencies/territories outside the 193 UN
  members (correctly out of scope), or ocean/Antarctica.
- **GADM fallback is non-commercial** — only reached if both UNOCHA and geoBoundaries lack a
  country/level; check `acquisition_manifest.csv` `source` column if license matters.

---

## 11. The agent-based SEIR model (`wwsim.abm`)

A planet-scale **agent-based** SEIR over the admin-2 nodes. Each (subsampled) person is one
agent with just `state` (uint8), `nodeid` (uint16), `timer` (uint8) = **4 bytes** — so full
resolution (~7.3 B agents ≈ 29 GB) fits in >32 GB RAM; `--subsample N` divides node
populations for light runs.

```shell
# 1/200-scale world, seed in China, intra-country spread only, 180 days:
python scripts/10_run_seir.py --subsample 200 --seed-iso3 CHN

# add the cross-border air network (top-250 airports):
python scripts/10_run_seir.py --subsample 200 --seed-iso3 CHN --air --air-top-n 250

# full resolution (needs >32 GB RAM and the top-N air file built in step 06):
python scripts/10_run_seir.py --subsample 1 --seed-iso3 CHN --nticks 365 --air --air-top-n 1000
```

| Flag | Meaning (default) |
|---|---|
| `--subsample` | divide each node's population (200) |
| `--nticks` | days (180) |
| `--beta` | transmission rate/day; `R0 = beta·infectious` (0.35) |
| `--incubation` / `--infectious` | E→I / I→R durations in days (4 / 6) |
| `--waning` | R→S days; 0 = plain SEIR (0) |
| `--intra-fraction` | fraction of a node's FOI exported within its country (0.1) |
| `--air` / `--air-top-n` / `--air-fraction` | enable cross-border air layer / airport cut-off / international export fraction |
| `--seed-iso3` / `--seed-nodeid` / `--seed-count` | where/how many initial infections |

**The two-tier coupling (the custom FOI).** Per tick, FOI is `beta · I/N` plus sparse
network coupling `+ Wᵀ·ft − ft·rowsum` (the sparse form of laser-generic's dense coupling):

- **Intra-country** `W` = each country's gravity matrix (block-diagonal) → within-country
  spread.
- **Air** `W` (optional) = the global cross-border network for the chosen top-N airports →
  only the gateway admin-2 nodes participate, and only across borders. Use the same `N` you
  built the air network with in step 06.

**Outputs** (`output/seir/`): `totals_<tag>.csv` (S/E/I/R + daily incidence),
`epicurve_<tag>.png`, `attack_<tag>.png` (per-admin-2 final attack-rate choropleth),
`countries_<tag>.png` (per-country arrival timing).

**Per-node history + animated choropleth.** Add `--save-history` to write the full per-node,
per-tick S/E/I/R to `output/seir/history_<tag>.npz`, then turn it into an mp4:

```shell
# 540-day run that feeds the animation:
python scripts/10_run_seir.py --subsample 20 --beta 0.35 --air --air-top-n 1000 \
    --seed-iso3 CHN --nticks 540 --tag world540_air --save-history

# 30-second movie of infectious prevalence (one frame per day; fps auto = ticks/seconds):
python scripts/12_animate_choropleth.py --history output/seir/history_world540_air.npz \
    --field I --seconds 30
```

`--field` selects `I`/`E`/`EI` (prevalence = compartment/population) or `attack` (cumulative
R/population). The script rasterizes the 45k polygons to a node-index grid once and recolors
it per frame, so a 540-frame movie renders in minutes. Output: `output/seir/anim_<tag>_<field>.mp4`.

**Validation.** A single well-mixed node reproduces the Kermack–McKendrick final size for
`R0 = beta·infectious_period` (`tests/test_abm.py`). Demo (1/200, seeded Shanghai): intra-
only stays confined to China (~14% global attack); adding the top-250 air layer crosses
borders and roughly doubles it (~29%).

**Notes / extending.** Durations are deterministic (clean R0); to randomize, sample
per-agent timers (e.g. via `laser.core.distributions`) when setting `timer`. The model is
closed (no births/deaths). The FOI lives entirely in `wwsim.abm.components.Transmission` —
the single place to change the epidemiology.
