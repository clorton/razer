# Estimate the agent capacity for a growing population reclaimed with [`squash()`](https://clorton.github.io/razer/reference/squash.md).

The mortality-aware companion to
[`calc_capacity()`](https://clorton.github.io/razer/reference/calc_capacity.md).
When dead agents' slots are reclaimed periodically with
[`squash()`](https://clorton.github.io/razer/reference/squash.md), the
slots needed are bounded by the **peak simultaneous living population**
— net births minus deaths — not by the cumulative number ever born
(which
[`calc_capacity()`](https://clorton.github.io/razer/reference/calc_capacity.md)
estimates for the no-reclaim case). This lets you model decades or
centuries without allocating one slot per agent ever born.

## Usage

``` r
calc_capacity_cdr(birthrates, deathrates, initial_pop, safety_factor = 1)
```

## Arguments

- birthrates:

  A 2-D `nsteps x nnodes` matrix of crude birth rates (per 1,000 per
  year), or a 2-D
  [Column](https://clorton.github.io/razer/reference/Column.md) (e.g.
  from
  [`values_map()`](https://clorton.github.io/razer/reference/values_map.md)).
  Each value in `[0, 100]`.

- deathrates:

  A 2-D `nsteps x nnodes` matrix (or 2-D
  [Column](https://clorton.github.io/razer/reference/Column.md)) of
  crude death rates, the same shape as `birthrates`. Each value in
  `[0, 100]`.

- initial_pop:

  A non-negative numeric vector of length `nnodes`: initial population.

- safety_factor:

  Non-negative headroom multiplier in `[0, 6]` (default `1`),
  controlling how much the death rate is underestimated. `0` credits
  deaths fully.

## Value

A numeric vector of length `nnodes` of estimated capacities
(whole-valued doubles). A `warning` is issued if any estimate exceeds
`.Machine$integer.max`.

## Details

The per-node daily birth and death rates
`lambda = (1 + rate/1000)^(1/365) - 1` are summed over all time steps;
the expected net-growth factor is `exp(sum lambda_b - sum lambda_d)`.
For a conservative bound the **death rate is underestimated** by the
safety factor — only a fraction `1 / (1 + safety_factor)` of the death
sum is credited, holding the rest back as headroom against a
lower-mortality (faster-growing) realization. So `safety_factor = 0`
credits deaths fully (the tightest, bare net-growth estimate); larger
values credit fewer deaths and reserve more slack. (Unlike
[`calc_capacity()`](https://clorton.github.io/razer/reference/calc_capacity.md),
no gross- births term enters — that would defeat the point of reclaiming
slots.)

## Errors

Stops if either rate grid is not 2-D, if their shapes or node counts
disagree with each other or with `length(initial_pop)`, if any
population is negative, if any rate is outside `[0, 100]`, or if
`safety_factor` is outside `[0, 6]`.

## See also

[`calc_capacity()`](https://clorton.github.io/razer/reference/calc_capacity.md)
(cumulative-births bound, no reclaim),
[`squash()`](https://clorton.github.io/razer/reference/squash.md).

## Examples

``` r
# one node, 100 years, CBR 30 / CDR 15 — peak-living bound for a squash-reclaimed run
br <- matrix(30, nrow = 100 * 365, ncol = 1)
dr <- matrix(15, nrow = 100 * 365, ncol = 1)
calc_capacity_cdr(br, dr, initial_pop = 1e6, safety_factor = 1)
#> [1] 9129894
```
