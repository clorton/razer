# Advance M→S, E→I, and I→`absorbing_state` for one tick — SIS / SIR / SEIS / SEIR.

`absorbing_state` is the untimed destination of the infectious period:
`S` (SIS/SEIS) or `R` (SIR/SEIR). Returns `list(waned, onset, cleared)`
of per-node counts (`cleared` is the I→`absorbing_state` flow).

## Usage

``` r
step_sir(state, timer, nodeid, count, n_nodes, inf_duration, absorbing_state)
```

## Arguments

- state, timer, nodeid, count, n_nodes:

  As in
  [`step_si()`](https://clorton.github.io/razer/reference/step_si.md).

- inf_duration:

  A Distribution for the infectious period set on E→I.

- absorbing_state:

  State code I clears to (`laser_states()[["S"]]` or `[["R"]]`).

## Value

`list(waned, onset, cleared)` of `integer[n_nodes]`.
