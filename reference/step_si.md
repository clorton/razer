# Advance M→S (maternal waning) and E→I (incubation) for one tick — SI / SEI.

`I` is terminal (no exit). Returns `list(waned, onset)` of per-node
counts.

## Usage

``` r
step_si(state, timer, nodeid, count, n_nodes, inf_duration)
```

## Arguments

- state:

  Per-agent `u8` state Column (mutated).

- timer:

  Per-agent `u16` timer Column (mutated; E→I draws an infectious timer).

- nodeid:

  Per-agent `u16` 0-based node-id Column.

- count:

  Number of active agents to process.

- n_nodes:

  Number of nodes (the length of each returned vector).

- inf_duration:

  A Distribution for the infectious period set on E→I.

## Value

`list(waned = integer[n_nodes], onset = integer[n_nodes])`.
