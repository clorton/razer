# Apply constant-population SIR vital dynamics for one tick.

`rate` is a per-node daily death-HAZARD-rate grid (a values map; the
caller converts a crude death rate of annual deaths per 1000 people to a
daily rate, e.g. `cdr / 1000 / 365`). For each of the first `count`
agents, this converts the node's rate to a probability `1 - exp(-rate)`
(once per node, like transmission) and, with that probability,
"unalives" the agent and replaces it with a newborn: the agent's `state`
is reset to Susceptible and its `timer` to 0.

## Usage

``` r
constant_pop_vitals_sir(
  state,
  timer,
  nodeid,
  count,
  rate,
  s_count,
  i_count,
  r_count,
  births,
  deaths,
  tick
)
```

## Arguments

- state:

  Per-agent `u8` state Column (mutated; deaths reset to Susceptible).

- timer:

  Per-agent `u16` countdown Column (mutated; deaths reset to 0).

- nodeid:

  Per-agent `u16` 0-based node-id Column.

- count:

  Number of active agents to process.

- rate:

  Per-node daily death-hazard-rate grid (`n_ticks x n_nodes`, from
  [`values_map()`](https://clorton.github.io/razer/reference/values_map.md));
  column `tick` is read.

- s_count, i_count, r_count:

  `n_ticks x n_nodes` i32 census Columns kept in sync (mutated at column
  `tick + 1`).

- births, deaths:

  `(n_ticks-1) x n_nodes` i32 flow Columns; column `tick` receives the
  per-node event counts (equal; mutated).

- tick:

  0-based tick index.

## Value

`NULL` (invisibly); the Columns are modified in place.

## Details

The S/I/R node census is updated IN PLACE at column `tick + 1` (the
working column the caller has already carried forward): a death out of I
decrements `I` and increments `S`; a death out of R decrements `R` and
increments `S`; a death out of S nets to zero. Every event (from any
state) is counted per node and written to BOTH the `births` and `deaths`
flow reports for `tick` (equal under constant population). Agents are
assumed to be in S, I, or R (it is the SIR variant). Parallelized with
private per-thread node buffers summed at the end.
