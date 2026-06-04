# Timer-based E→I transition (SEIR kernel).

For each exposed agent, decrements `timer` by 1. When `timer` reaches 0
the agent transitions to state I and `timer` is set to a fresh draw from
`inf_dist`, the infectious-period distribution. Pass `dist_constant(d)`
for a fixed period of `d` ticks, or e.g. `dist_normal(mean, variance)`
for a stochastic per-agent period. Draws are rounded to the nearest tick
and clamped to a minimum of 1.

## Usage

``` r
step_exposed_ei(people, inf_dist)
```

## Arguments

- people:

  LaserFrame of agents.

- inf_dist:

  A `Distribution` (e.g. from
  [`dist_constant()`](https://clorton.github.io/razer/reference/dist_constant.md)
  or
  [`dist_normal()`](https://clorton.github.io/razer/reference/dist_normal.md))
  giving the infectious period in ticks; sampled and written to `timer`
  on E→I.

## Details

A fixed-state shorthand for
[`step_timer_expire_set()`](https://clorton.github.io/razer/reference/step_timer_expire_set.md)`(people, E, I, inf_dist)`.

**RNG:** thread-local — each Rayon worker draws from its own
`thread_rng` (Pattern B: the kernel owns the RNG and passes it into the
sampler). The single `inf_dist` handle is shared across threads by
reference.
