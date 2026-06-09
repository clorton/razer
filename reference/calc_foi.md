# Compute the per-node force of infection (FOI) for one tick.

Computes the frequency-dependent, network-redistributed FOI **rate**
into column `tick` of `foi`. The local rate per node is
`r[k] = beta[k] * seasonality[k] * infected[k] / population[k]`,
redistributed as
`foi[k] = r[k] * (1 - sum_j W[k, j]) + sum_i r[i] * W[i, k]`. The
transmission kernels turn this rate into a per-tick probability
`1 - exp(-foi)`.

## Usage

``` r
calc_foi(infected, population, beta, seasonality, network, foi, tick)
```

## Arguments

- infected:

  Infectious-count census Column (`nodes$I`), `slice_len == n_nodes`.

- population:

  Per-node population census Column (`nodes$N`); the FOI denominator.

- beta:

  Transmission-coefficient grid (`n_ticks x n_nodes`, from
  [`values_map()`](https://clorton.github.io/razer/reference/values_map.md)).

- seasonality:

  Seasonal-modifier grid (`n_ticks x n_nodes`, from
  [`values_map()`](https://clorton.github.io/razer/reference/values_map.md)).

- network:

  The `n_nodes x n_nodes` coupling weights, column-major. Either an R
  numeric matrix, OR a razer
  [Column](https://clorton.github.io/razer/reference/Column.md) of
  `n_nodes * n_nodes` f64 holding the matrix in column-major order. The
  Column form avoids re-marshalling the matrix from R on every tick
  (build it once, e.g.
  `nc <- allocate_vector("f64", n, n); nc$set(as.vector(W))`).

- foi:

  A 2-D f64 Column (`(n_ticks-1) x n_nodes`); column `tick` is
  overwritten.

- tick:

  0-based tick index: reads `beta[tick]`/`seasonality[tick]` and
  `infected[tick]`/`population[tick]` (the start-of-interval census),
  writes `foi[tick]`.

## Value

`NULL` (invisibly); the result is written into `foi`.

## Details

Index conventions: `infected` and `population` are census buffers read
at the **start-of-interval** column `tick` (the settled census recorded
for tick `tick`, the infectious count and population at the start of the
interval `tick → tick+1`); `beta` and `seasonality` are exogenous
modifier grids read at the interval column `tick`; the result is written
to `foi[tick]`.

**Ordering and the effective infectious period.** Because `calc_foi`
reads the SETTLED census at column `tick` — not the working column
`tick+1` that this interval's step kernel and transmission build — its
result does NOT depend on where the step kernel runs. That lets every
model use ONE ordering:
`carry_forward → step → calc_foi → transmission`, with `calc_foi` placed
immediately before `transmission`. Under it an agent contributes to the
FOI on exactly the `D` census columns it occupies (it enters `I` at
column `entry+1` via either the step kernel's E→I or transmission's S→I,
and recovers `D` columns later), so the realized basic reproduction
number is the full `R0 = beta * D` — for both direct S→I (SIR) and
SEIR-style entry, with no `beta*(D-1)` artifact and no per-family
special-casing.
