# Pipeline & data model

## Stages

```
01 acquire_countries ─┐                         laser-init + rasterio fallback
   (per country)       ├─► data/countries/<ISO3>/<ISO3>_admin2.gpkg
01b fix_failed ───────┘                         + output/nodes/acquisition_manifest.csv

02 build_nodes  ──► output/nodes/global_admin2_nodes.{gpkg,parquet}   (the node table)

03 ingest_flights ─► data/interim/airports_master.parquet
                     data/interim/airport_edges.parquet               (airport O-D + proxy)

04 assign_airports ─► data/interim/airport_node_assignment.parquet    (airport → admin-2 node)

05 intra_country_matrices ─► output/networks/<ISO3>_gravity.{npz,parquet}  (per-country)

06 global_air_network ─────► output/networks/global_air_network.{parquet,npz}  (cross-border)

07 combine_network ────────► output/networks/global_combined_network.npz   (intra + inter)

08 plots ──────────────────► output/plots/*.png
```

## The node table (the spine)

`output/nodes/global_admin2_nodes.parquet` — one row per admin-2 unit worldwide:

| column | meaning |
|---|---|
| `global_nodeid` | dense 0..N-1 id; **the index for every matrix** |
| `iso3` | country (UN member) |
| `adm2_name`, `adm2_id` | admin-2 label and stable source id (pcode/shapeID/GID_2) |
| `adm1_name` | parent admin-1 where the source provides it |
| `population` | 2015 1km WorldPop, summed within the polygon |
| `lon`, `lat` | polygon centroid (gravity distances) |
| `geometry` | admin-2 polygon (EPSG:4326) |

Every network indexes rows/cols by `global_nodeid`, so all matrices share one coordinate
system of nodes.

## The networks

- **Per-country (intra) gravity** — `w_ij = k·P_i^a·P_j^b / D_ij^c` over each country's
  admin-2 nodes. One matrix per country; together they form the **block-diagonal**
  intra-country migration.
- **Global (inter) air** — directed admin-2↔admin-2 edges built **only between countries**,
  from the **top-N airports** by passenger volume, with **all airports in one admin-2
  aggregated** into that node. Weight = summed seat-capacity proxy.
- **Combined** — `gravity_scale · intra ⊕ air_scale · inter` (+ optional rail), one global
  sparse matrix; the inter-nation coupling for a worldwide metapopulation SEIR run.

## Using it in LASER / Razer

Load the node table for per-patch population and coordinates, and the combined matrix
(`scipy.sparse`) as the migration/mixing network. Node ids align across the node table and
all `.npz` matrices:

```python
import geopandas as gpd, scipy.sparse as sp
nodes = gpd.read_parquet("output/nodes/global_admin2_nodes.parquet")
M = sp.load_npz("output/networks/global_combined_network.npz")   # (N, N), N == len(nodes)
pop = nodes["population"].to_numpy()
```

`M[i, j]` is the directed coupling from node `i` to node `j` (intra-country gravity on the
diagonal blocks, cross-border air off the blocks). Row-normalize for migration
probabilities, or scale into a force-of-infection mixing term, per your model's convention.
