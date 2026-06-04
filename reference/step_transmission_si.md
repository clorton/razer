# Stochastic S→I transmission step (SI / SIR kernel).

**Parallelism:** the active agent array is split into fixed chunks — one
per Rayon worker thread — matching the `prange` pattern. The aggregation
phase uses a thread-local fold then reduces by element-wise summation.
The FOI phase runs each chunk independently; RNG is thread-local
(`rand::thread_rng()`).

## Usage

``` r
step_transmission_si(people, nodes, beta, inf_dist)
```

## Arguments

- people:

  LaserFrame of agents.

- nodes:

  LaserFrame of patches/nodes.

- beta:

  Transmission rate (force of infection per infectious contact per
  tick).

- inf_dist:

  A `Distribution` (e.g.
  [`dist_constant()`](https://clorton.github.io/razer/reference/dist_constant.md)
  or
  [`dist_normal()`](https://clorton.github.io/razer/reference/dist_normal.md))
  giving the infectious period in ticks; sampled per newly infected
  agent and written to `timer` on S→I (rounded to whole ticks, clamped
  to a minimum of 1).

## Details

All operations are performed in-place on the backing arrays; no copies
are made.

Node-level I counts are computed from current agent states and written
to `nodes$I` (overwriting the previous value).

**Required people properties:** `state` (integer), `node` (integer,
0-based), `timer` (integer). **Required nodes properties:** `N`
(integer, total population per node), `I` (integer, will be
overwritten).
