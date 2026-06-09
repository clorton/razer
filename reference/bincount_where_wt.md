# Sum a weight, per group, over the agents whose property satisfies a comparison.

The weighted twin of
[`bincount_where()`](https://clorton.github.io/razer/reference/bincount_where.md):
for each group `g` in `0..n_groups`, sums `weights[i]` over the first
`count` agents that both have `group[i] == g` AND satisfy
`prop[i] <op> value`. Use it for predicate-filtered weighted aggregates
by node — e.g. total infectiousness of symptomatic agents per node, or
person-days under five — in one parallel pass with no copy of `prop` or
`weights` into R.

## Usage

``` r
bincount_where_wt(
  group,
  n_groups,
  prop,
  op,
  value,
  weights,
  count,
  counts = NULL,
  slot = 0L
)
```

## Arguments

- group:

  An integer-typed
  [Column](https://clorton.github.io/razer/reference/Column.md) of group
  indices (`i8`..`u32`) — e.g. `nodeid`.

- n_groups:

  Number of groups; a non-negative integer.

- prop:

  A numeric
  [Column](https://clorton.github.io/razer/reference/Column.md) holding
  the per-agent property to test.

- op:

  Comparison string: one of `"eq"`, `"ne"`, `"lt"`, `"le"`, `"gt"`,
  `"ge"`.

- value:

  The threshold the property is compared against.

- weights:

  A numeric
  [Column](https://clorton.github.io/razer/reference/Column.md) of
  per-agent weights to sum (any numeric type).

- count:

  How many leading agents to scan — the ACTIVE agent count, typically
  `people$count`. Required (see
  [`bincount_where()`](https://clorton.github.io/razer/reference/bincount_where.md)):
  scanning the Column's full length would tally reserved, inactive slots
  in an over-allocated population.

- counts:

  Optional numeric
  [Column](https://clorton.github.io/razer/reference/Column.md) to
  receive the sums; when omitted a numeric result vector is allocated
  and returned instead.

- slot:

  Which slice of `counts` to write when `counts` is supplied. Defaults
  to `0`.

## Value

When `counts` is `NULL`, a numeric vector of per-group sums; otherwise
`NULL` invisibly (the result is written into `counts`).

## Details

Two output modes, as in
[`bincount_where()`](https://clorton.github.io/razer/reference/bincount_where.md):
leave `counts` `NULL` (the default) and a numeric vector of length
`n_groups` is allocated and returned; or pass a numeric
[Column](https://clorton.github.io/razer/reference/Column.md) `counts`
and the per-group sums are written into its slice `slot`.

## See also

[`bincount()`](https://clorton.github.io/razer/reference/bincount.md),
[`bincount_wt()`](https://clorton.github.io/razer/reference/bincount_wt.md),
[`bincount_where()`](https://clorton.github.io/razer/reference/bincount_where.md).

## Examples

``` r
states  <- laser_states()
state   <- allocate_scalar("u8",  5L); state$set(c(2, 2, 0, 2, 1))     # I I S I E
#> NULL
nodeid  <- allocate_scalar("u16", 5L); nodeid$set(c(0, 0, 1, 1, 1))
#> NULL
shed    <- allocate_scalar("f64", 5L); shed$set(c(1.0, 0.5, 9, 2.0, 9))# infectiousness
#> NULL
# total shedding of infectious (state==I) agents per node:
bincount_where_wt(nodeid, 2L, state, "eq", states[["I"]], shed, count = 5L)  # 1.5  2.0
#> [1] 1.5 2.0
```
