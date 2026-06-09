# Changelog

All notable changes to the `worldwide` (`wwsim`) data pipeline are recorded here.

## [Unreleased]

### Added — technical write-up
- `writeup.md`: a thorough technical report (data collection, processing, the agent-based
  SEIR engine, laptop-feasibility choices, the Numba kernels, gravity↔air interplay, tests,
  performance) with mermaid diagrams and figures. `scripts/15_writeup_media.py` regenerates
  the figures into `media/` (new analytical SVGs: gravity decay, memory budget, admin-2
  distribution, single-node SEIR validation, global epicurve; plus copied choropleths and an
  animation still). Note: the COVID-comparison material is intentionally excluded from the
  write-up.

### Fixed
- `WorldSEIR` node-id initialization no longer builds a multi-GB `int64` temporary
  (`np.repeat(np.arange(..., int64), pops)` was ~29 GB at subsample 2 / ~3.66 B agents,
  ~3x the steady-state, OOM-killing in containers with a hard cgroup limit and little swap —
  macOS hid it via memory compression + swap). Node ids are now filled block-by-block into
  the preallocated `uint16` column (no large temporary); init peak drops ~50 GB → ~16 GB.

### Added — per-country Rt forcing (real COVID transmissibility into the model)
- `covid_reference.download_owid_full` / `country_rt_matrix` / `rt_forcing`: build a
  per-country, per-tick transmission multiplier `m_c(t) = Rt_c(t) / baseline` from OWID's
  `reproduction_rate` (Arroyo-Marioli Kalman estimate; 191/193 UN members, 2020-01-23 →
  2023-01-02). Sim day 0 maps to a calendar anchor; missing Rt → `m=1` (no forcing).
- `WorldSEIR(..., forcing=(m, iso3_order))` + `Transmission`: the FOI uses
  `beta · m_country(node)(t) · I/N`, so each admin-2 node applies its country's factor (and
  the multiplier propagates into the network coupling via `ft0`).
- `scripts/10_run_seir.py --rt-forcing [--rt-anchor --rt-baseline --rt-clip-max]`. Tests in
  `tests/test_abm.py` (unity-forcing is a no-op; zero halts transmission; >1 raises attack).

### Added
- Project scaffolding: `wwsim` package (`src/wwsim`), `pyproject.toml`, `.venv` (uv, Python 3.12).
- `wwsim.logging`: centralized package logger (`from wwsim.logging import logger`).
- `wwsim.config`: `Config`/`GravityParams`/`AirNetworkParams` dataclasses, on-disk layout,
  and `load_config()` YAML overlay.
- `wwsim.countries`: canonical list of the 193 UN member states (ISO3) with guards.
- Vendored `laser-init` (from github.com/laser-base/laser-init) installed editable into
  `.venv`; validated end-to-end on Comoros (`COM_admin2.gpkg`, 17 admin-2 units, 1km
  WorldPop 2015 constrained UN-adjusted, EPSG:4326).
- `wwsim.acquire`: `Acquirer` driving laser-init per country with a source waterfall
  (unocha→geoboundaries→gadm) and admin-2→admin-1 fallback; UNOCHA global geodatabase is
  read once and cached in memory (not re-read per country). Resumable; writes a manifest.
- `wwsim.zonal` + `Acquirer.acquire_country_rasterio`: rasterio zonal-sum fallback for
  countries whose WorldPop file is a global canvas (x0=-180) that RasterToolkit rejects
  (RUS/FJI/KIR/TUV antimeridian + ERI).
- `wwsim.nodes`: fuse per-country gpkgs into one global admin-2 node table
  (`global_nodeid, iso3, adm2_name, adm2_id, adm1_name, population, lon, lat, geometry`);
  harmonizes UNOCHA/geoBoundaries/GADM schemas; repairs name mojibake.
- `wwsim.flights`: OurAirports + OpenFlights ingestion; curated aircraft seat table
  (≈97% equipment-code coverage); directed airport edges with a seat-capacity passenger
  proxy.
- `wwsim.airports_assign`: airport → admin-2 node by point-in-polygon, with border
  disambiguation and same-country nearest-node fallback.
- `wwsim.gravity`: per-country gravity migration matrices over admin-2 nodes (haversine
  distances), `.npz` + Parquet outputs, block-diagonal global assembly.
- `wwsim.air_network`: one global cross-border air network over admin-2 nodes — top-N
  airports, cross-border-only, airports aggregated per admin-2.
- `wwsim.network` + `wwsim.rail` + `wwsim.oag`: pluggable multi-modal combiner
  (`ModeNetwork`/`combine_modes`); rail adapter + documented sources; OAG drop-in adapter.
- `wwsim.plotting`: global population choropleth, air-network map, airport-assignment map,
  edge-weight distribution, per-country gravity map.
- `scripts/01..08`, `01b_fix_failed`, `run_all`: end-to-end pipeline runners.
- Test suite (`tests/`, given-when-then): 33 tests across countries, config, gravity,
  flights, airports-assign, air-network, network, nodes, rail/oag — all passing.
- `docs/datasources.md`, `docs/pipeline.md`: sources/licenses and the data model.
- `--top-n` flag on steps 06 and 07: build the air network and the combined matrix for one
  or more airport cut-offs; multiple values write suffixed `*_top<N>.{parquet,npz}` files so
  cut-offs coexist (`0` = no cap). `save_air_network`/`save_combined_network` take a `stem`.
- `README.md` + `USAGE.md`: full step-by-step usage guide.

### Added — agent-based SEIR model (`wwsim.abm`)
- `wwsim.abm.model.WorldSEIR` / `SEIRParams`: agent-based SEIR over all admin-2 nodes on
  `laser-core` `LaserFrame`s. Parsimonious **4 bytes/agent** — `state` uint8, `nodeid`
  uint16, `timer` uint8 (one timer reused for incubation/infectious/waning) — so full
  resolution (~7.3B agents ≈ 29 GB) fits in >32 GB; `subsample` divides node populations
  for light test runs.
- `wwsim.abm.components`: single-pass **Progression** (branch-on-entry-state E→I→R, optional
  R→S waning) and a **custom Transmission/FOI** that replaces laser-generic's dense
  `ft[:,None]*network` with **sparse** matvecs (`ft + Wᵀ·ft − ft·rowsum`). Razer ordering
  (progression before transmission) yields exact `R0 = beta · infectious_period`.
- `wwsim.abm.networks`: sparse, row-normalized **intra-country (block-diagonal)** coupling
  and optional **cross-border air** coupling (only gateway admin-2 nodes of the chosen
  top-N airports participate).
- `wwsim.abm.plots`, `scripts/10_run_seir.py`: run + epicurve / attack-rate choropleth /
  per-country curves. Tests `tests/test_abm.py` (6, incl. Kermack–McKendrick validation).
- **Validation**: single well-mixed node attack fraction 0.737 vs KM 0.732 (R0=1.8).
- **Demo** (1/200 scale, seeded Shanghai, 180 d, β=0.35): intra-only confines to the seed
  country (final attack **14.4%**); adding the air network (top-250) crosses borders via 237
  gateway nodes → final attack **28.6%**, later/larger peak. ~13 s/run (45,406 nodes, 36.6M
  agents).

### Added — per-node history + animated choropleth
- `WorldSEIR.save_history()` + `10_run_seir.py --save-history`: write per-node, per-tick
  S/E/I/R to `output/seir/history_<tag>.npz` (compressed; ~100 MB for 45,406 nodes × 541
  ticks).
- `scripts/12_animate_choropleth.py`: render an mp4 choropleth movie from a history file —
  rasterizes the admin-2 polygons to a node-index grid **once** and recolors per frame
  (`--field I|E|EI|attack`, `--seconds`, `--fps`, `--width`). 540-frame movie renders in
  minutes, not hours.
- 540-day demo (subsample 20, β=0.35, air top-1000, seeded Shanghai): 365.8M agents, ran in
  ~2 min; peak I 19.2M (day 221), final attack 84.5%; 30-s animation at 18 fps.

### Added — real COVID reference data (`wwsim.covid_reference`)
- `wwsim.covid_reference` + `scripts/11_covid_reference.py`: download Our World in Data
  cases/deaths (2020-01-05 → 2024-08-04; COVID-19 did not exist in 2000-2004) and plot the
  global curve (daily + cumulative cases, daily deaths), per-country curves, and a
  normalized single-wave shape overlay vs a `wwsim` SEIR run. Outputs in `output/covid/`.
  Real totals: 775.9M reported cases, 7.06M deaths.
- `covid_reference.country_value_matrix` + `owid_to_iso3` (193/193 UN members mapped) and
  `scripts/13_animate_covid.py`: a **country choropleth animation** of real COVID (daily new
  cases/deaths per 100k, 7-day avg). Reuses the admin-2 node raster colored by each node's
  country value (no polygon dissolve). Default 730 days (2020-01-05 → 2022-01-03) in 30 s.
- `wwsim.choropleth_anim` (shared renderer: `rasterize_nodes`, `country_vmax`,
  `render_country_movie`) + `scripts/14_animate_sim_cases.py`: render the **simulation as
  daily-new-infections per 100k**, country-aggregated, on the **identical color scale, norm,
  colormap, and aggregation** as the COVID movie (shared `vmax` = the COVID series' p99), so
  the two mp4s are directly comparable. Sim incidence is derived from the saved history as
  `S[t]-S[t+1]` (exact in a no-waning SEIR).

### Build results (first full run, 2015)
- **193/193 UN member states** acquired (105 via UNOCHA COD, 79 via geoBoundaries, plus
  GADM/fallbacks; 174 at admin-2, 14 at admin-1). 5 antimeridian/global-canvas countries
  (RUS, FJI, KIR, TUV, ERI) repaired via the rasterio fallback; KIR via the 2000-2020
  raster product.
- **45,406 admin-2 nodes**, total population **7.32 billion** (~99.5% of the 2015 world).
- Air: 9,056 airports, 67,663 routes → **15,221 cross-border admin-2 edges** (top-1000
  airports). Top corridors GBR↔ESP, DEU↔ESP, USA↔GBR, MEX↔USA (match real busiest routes).
- Combined global migration matrix **45,406 × 45,406**, 72.8M nonzeros (intra gravity
  block-diagonal + inter air). Total on-disk ~8.4 GB (budget 150 GB).
