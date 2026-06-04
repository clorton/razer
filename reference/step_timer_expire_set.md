# Generalized timer-expiry transition into a state that has *its own* duration.

For each agent currently in `from_state`, decrements `timer` by 1; when
the timer reaches 0 the agent moves to `to_state` and `timer` is reset
to a fresh per-agent draw from `duration_dist` (rounded to whole ticks,
clamped to a minimum of 1). This is the "transition to a state with its
own duration timer" generalization (laser-generic's
`nb_timer_update_timer_set`), e.g. E→I in SEIR or I→R with waning in
SEIRS.

## Usage

``` r
step_timer_expire_set(people, from_state, to_state, duration_dist)
```

## Arguments

- people:

  LaserFrame of agents.

- from_state:

  Integer state code an agent must currently occupy to be eligible.

- to_state:

  Integer state code an agent moves to when its timer expires.

- duration_dist:

  A `Distribution` (e.g.
  [`dist_constant()`](https://clorton.github.io/razer/reference/dist_constant.md)
  or
  [`dist_normal()`](https://clorton.github.io/razer/reference/dist_normal.md))
  giving the destination state's duration in ticks; sampled per
  transitioning agent and written to `timer`.

## Details

This is the engine behind
[`step_exposed_ei()`](https://clorton.github.io/razer/reference/step_exposed_ei.md)
(E→I) and
[`step_infectious_ir()`](https://clorton.github.io/razer/reference/step_infectious_ir.md)
(I→R):
`step_timer_expire_set(people, laser_states()[["E"]], laser_states()[["I"]], inf_dist)`
is exactly `step_exposed_ei(people, inf_dist)`.

**RNG:** thread-local (Pattern B) — each Rayon worker draws from its own
`thread_rng`; the single `duration_dist` handle is shared across threads
by reference.

**Required people properties:** `state`, `timer` (both integer).
