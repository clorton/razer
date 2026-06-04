# Timer-based I→R transition (SIR / SEIR / SEIRS kernel).

For each infectious agent, decrements `timer` by 1. When `timer` reaches
0 the agent transitions to state R and `timer` is set to a draw from
`imm_dist`, the immunity (waning) period.

## Usage

``` r
step_infectious_ir(people, imm_dist)
```

## Arguments

- people:

  LaserFrame of agents.

- imm_dist:

  A `Distribution` (e.g.
  [`dist_constant()`](https://clorton.github.io/razer/reference/dist_constant.md)
  or
  [`dist_normal()`](https://clorton.github.io/razer/reference/dist_normal.md))
  giving the immunity period in ticks; sampled and written to `timer` on
  I→R (rounded to whole ticks, clamped to a minimum of 1). Use
  `dist_constant(0)` for SIR / SEIR.

## Details

A fixed-state shorthand for
[`step_timer_expire_set()`](https://clorton.github.io/razer/reference/step_timer_expire_set.md)`(people, I, R, imm_dist)`.

For SEIRS, pass the desired immunity distribution so that
`step_recovered_rs` counts down R→S from a fresh per-agent draw. For SIR
and SEIR (no waning, `step_recovered_rs` never called) the R timer is
never read, so any distribution works — `dist_constant(0)` is the
conventional "don't care".

**RNG:** thread-local — each Rayon worker draws from its own
`thread_rng` (Pattern B). The single `imm_dist` handle is shared across
threads by reference.
