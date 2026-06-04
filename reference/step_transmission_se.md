# Stochastic S→E exposure step (SEIR kernel).

Same FOI computation and parallelism as `step_transmission_si`, but
newly exposed agents move to state E and `timer` is set to a draw from
`exp_dist` (incubation period). Pair with `step_exposed_ei` to complete
E→I.

## Usage

``` r
step_transmission_se(people, nodes, beta, exp_dist)
```

## Arguments

- people:

  LaserFrame of agents.

- nodes:

  LaserFrame of patches/nodes.

- beta:

  Transmission rate.

- exp_dist:

  A `Distribution` (e.g.
  [`dist_constant()`](https://clorton.github.io/razer/reference/dist_constant.md)
  or
  [`dist_normal()`](https://clorton.github.io/razer/reference/dist_normal.md))
  giving the incubation period in ticks; sampled per newly exposed agent
  and written to `timer` on S→E (rounded to whole ticks, clamped to a
  minimum of 1).

## Details

**Required people properties:** `state`, `node`, `timer` (all integer).
**Required nodes properties:** `N`, `I` (integer; `I` is overwritten).
