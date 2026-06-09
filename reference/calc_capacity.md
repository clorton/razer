# Estimate the agent capacity needed for a growing population.

Projects each node's population forward under its (possibly
time-varying) crude birth rate and returns a per-node capacity to
preallocate, inflated by an optional safety factor. The births are
treated as geometric growth: a daily growth rate
`lambda = (1 + CBR/1000)^(1/365) - 1` is summed over all time steps and
exponentiated (`exp(sum lambda)`), a geometric-Brownian-motion-style
estimate of expected growth. The safety factor adds a multiple of the
growth's standard-deviation-like term `sqrt(exp(sum lambda)) - 1` as
headroom for stochastic variation.

## Usage

``` r
calc_capacity(birthrates, initial_pop, safety_factor = 1)
```

## Arguments

- birthrates:

  A 2-D `nsteps x nnodes` numeric matrix of crude birth rates (births
  per 1,000 individuals per year), or a 2-D razer
  [Column](https://clorton.github.io/razer/reference/Column.md) (e.g.
  from
  [`values_map()`](https://clorton.github.io/razer/reference/values_map.md))
  whose `$values()` is such a matrix. Each value must be in `[0, 100]`.

- initial_pop:

  A numeric vector of length `nnodes`: the initial population per node.
  Must be non-negative.

- safety_factor:

  Non-negative headroom multiplier in `[0, 6]` (default `1`). `0` gives
  the bare expected-growth estimate; larger values reserve more slack.

## Value

A numeric vector of length `nnodes` of estimated capacities
(whole-valued doubles, which represent integers exactly up to `2^53`). A
`warning` is issued if any estimate exceeds `.Machine$integer.max` (R's
32-bit signed integer max) — the largest count razer's allocators
(`allocate_scalar`/`allocate_vector`, whose `count` is an `i32`) can
accept.

## Errors

Stops if `birthrates` is not 2-D, if its node count does not match
`length(initial_pop)`, if any population is negative, if any birth rate
is outside `[0, 100]`, or if `safety_factor` is outside `[0, 6]`.

## Examples

``` r
# two nodes, one year of a constant CBR of 40 per 1,000
br <- matrix(40, nrow = 365, ncol = 2)
calc_capacity(br, initial_pop = c(1e6, 5e5), safety_factor = 1)
#> [1] 1060599  530300
```
