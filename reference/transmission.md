# Stochastic transmission S→`to_state`, returning new infections per node.

Converts column `tick` of `foi` to a per-node probability
`1 - exp(-foi)` (once per node), then for each susceptible agent moves
it — with that probability — into `to_state` (`E` or `I`), setting its
u16 `timer` from `duration` (the incubation or infectious period). The
node census is NOT touched: the per-node count of new infections is
RETURNED, and the caller applies the `S` ↓ / `to_state` ↑ delta and
records incidence as its model requires (it knows whether this is S→E or
S→I).

## Usage

``` r
transmission(state, timer, nodeid, count, foi, tick, to_state, duration)
```

## Arguments

- state:

  Per-agent `u8` state Column (mutated).

- timer:

  Per-agent `u16` timer Column (mutated; set from `duration`).

- nodeid:

  Per-agent `u16` 0-based node-id Column.

- count:

  Number of active agents to process.

- foi:

  `n_ticks x n_nodes` f64 FOI Column (from
  [`calc_foi()`](https://clorton.github.io/razer/reference/calc_foi.md));
  column `tick` read.

- tick:

  0-based tick index.

- to_state:

  State code new infections enter (`laser_states()[["E"]]` or
  `[["I"]]`).

- duration:

  A Distribution from which the receiving state's timer is drawn.

## Value

An integer vector of new infections per node (length `n_nodes`).
