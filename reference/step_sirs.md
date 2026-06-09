# Advance M→S, E→I, I→R (with waning immunity), and R→S for one tick — SIRS / SEIRS.

On I→R a fresh immunity timer is drawn from `imm_duration`; R→S fires
when it expires. Returns `list(waned, onset, recovered, waned_r)` of
per-node counts.

## Usage

``` r
step_sirs(state, timer, nodeid, count, n_nodes, inf_duration, imm_duration)
```

## Arguments

- state, timer, nodeid, count, n_nodes:

  As in
  [`step_si()`](https://clorton.github.io/razer/reference/step_si.md).

- inf_duration:

  A Distribution for the infectious period set on E→I.

- imm_duration:

  A Distribution for the immunity period set on I→R.

## Value

`list(waned, onset, recovered, waned_r)` of `integer[n_nodes]`.
