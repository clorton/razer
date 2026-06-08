# Data sources & licenses

All inputs are open and target **vintage ~2015**. Provenance for every downloaded file is
recorded by `laser-init` in `~/.laser/cache/provenance.json` and per-output
`provenance.json`.

## Administrative boundaries (admin-2)

Acquired via [`laser-init`](https://github.com/laser-base/laser-init) with a per-country
**source waterfall** (`wwsim.acquire`): UNOCHA first (authoritative where it reaches
admin-2), then geoBoundaries (global coverage), then GADM; and **admin-2 → admin-1** level
fallback for the ~dozen microstates without admin-2 anywhere.

| Source | What | Coverage at ADM2 | License |
|---|---|---|---|
| **UNOCHA COD-AB** | One global "matched" geodatabase (`global_admin_boundaries_matched_latest.gdb`), filtered per country. Fields `adm{0..2}_name`, `adm2_pcode`, `iso3`. | ~humanitarian-priority subset (≈90 countries with ADM2 in the matched gdb) | Per-country (mostly open government); verify per dataset |
| **geoBoundaries gbOpen v6.0.0** | Per-country ADM2 shapefile zips. Fields `shapeName`, `shapeID`, `shapeGroup`. | ~180 of 193 | CC-BY 4.0 (attribution required) |
| **GADM 4.1** | Per-country shapefile zips. Fields `NAME_2`, `GID_2`. | global | **Academic / non-commercial only** — used only as a last-resort fallback |

> The actual source used per country is recorded in `output/nodes/acquisition_manifest.csv`
> (`source`, `level` columns).

## Population raster

| Source | Product | Resolution / year | License |
|---|---|---|---|
| **WorldPop** | Global 2015–2030 R2025A, **constrained, UN-adjusted** (`{iso}_pop_2015_CN_1km_R2025A_UA_v1.tif`) | **1 km**, **2015**, persons/pixel, EPSG:4326 | CC-BY 4.0 |

Population is aggregated to admin-2 polygons by **RasterToolkit** (`raster_clip`, MIT). For a
handful of countries whose raster is delivered on a full global canvas (tie point
`x0 = -180`, which trips RasterToolkit's strict bounds assertion) — the antimeridian
crossers **RUS, FJI, KIR, TUV** and the quirk **ERI** — `wwsim.zonal` sums population with
**rasterio** instead (windowed per polygon). Same numbers, no assertion.

## Air travel (open OAG substitute)

No free dataset replaces OAG passenger origin-destination *volumes*, so we approximate:

| Source | What | License |
|---|---|---|
| **OurAirports** `airports.csv` | Airport master: coordinates, ISO country, IATA/ICAO. | Public domain |
| **OpenFlights** `routes.dat` (~2014) | Carrier-level origin-destination routes + aircraft `equipment` codes. | ODbL (attribute OpenFlights; share-alike on redistributed derived DB) |
| `wwsim.flights.AIRCRAFT_SEATS` | Curated aircraft-type → typical-seats table (covers ~97% of route equipment codes). | This repo |
| **World Bank** `IS.AIR.PSGR` | Country-level passengers carried (calibration only; no O-D). | CC-BY 4.0 |

**Passenger-volume proxy.** Each carrier-route is weighted by representative **seat
capacity** (mean seats of its aircraft types); summing across carriers on a directed
airport pair gives a relative "seats offered" volume. Alternative: route multiplicity
(`passenger_proxy="routes"`). A licensed OAG O-D table drops in via `wwsim.oag` with no
downstream change.

## Demographics (optional, produced by laser-init)

| Source | What | License |
|---|---|---|
| **UN World Population Prospects 2024** | CBR/CDR, age distribution, life tables (per country, 2015). | CC-BY 3.0 IGO |

These are written per country by the full `laser-init` CLI but are **not required** by the
network pipeline; `wwsim.acquire` skips them for speed.

## Rail (documented, pluggable — not built by default)

See `wwsim.rail`. Candidate open sources: OpenStreetMap / OpenRailwayMap (ODbL), Eurostat
`rail_pa_*` (EU cross-border flows), UIC statistics, national GTFS. Any of these can be
reduced to an admin-2 edge list and summed into the global network as another
`ModeNetwork`.

## Attribution string (for derived outputs)

> Boundaries: UNOCHA COD-AB / geoBoundaries (CC-BY 4.0) / GADM (non-commercial).
> Population: WorldPop (CC-BY 4.0). Aggregation: RasterToolkit (MIT).
> Air routes: OpenFlights (ODbL); airports: OurAirports (public domain).
> Demographics: UN WPP 2024 (CC-BY 3.0 IGO).
