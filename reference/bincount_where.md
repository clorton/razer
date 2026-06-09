# Count, per group, the agents whose property satisfies a comparison.

A predicate-filtered, count-aware
[`bincount()`](https://clorton.github.io/razer/reference/bincount.md):
for each group `g` in `0..n_groups`, counts how many of the first
`count` agents both have `group[i] == g` AND satisfy
`prop[i] <op> value`. This answers flexible agent queries directly on
the Columns — e.g. "exposed by node" (`prop = state`, `op = "eq"`,
`value = laser_states()[["E"]]`) or "under-fives by node" (`prop = dob`,
`op = "gt"`, `value = tick - 5 * 365`, since `dob` is the negative age)
— in one parallel pass with no copy of `prop` into R.

## Usage

``` r
bincount_where(
  group,
  n_groups,
  prop,
  op,
  value,
  count,
  counts = NULL,
  slot = 0L
)
```

## Arguments

- group:

  An integer-typed
  [Column](https://clorton.github.io/razer/reference/Column.md) of group
  indices (`i8`..`u32`), each in `0..n_groups` — typically `nodeid`.

- n_groups:

  Number of groups; a non-negative integer.

- prop:

  A numeric
  [Column](https://clorton.github.io/razer/reference/Column.md) holding
  the per-agent property to test (compared as a double).

- op:

  Comparison string: one of `"eq"`, `"ne"`, `"lt"`, `"le"`, `"gt"`,
  `"ge"`.

- value:

  The threshold the property is compared against.

- count:

  How many leading agents to scan — the ACTIVE agent count, typically
  `people$count`. Required: there is no implicit population size, and
  the Column may be over-allocated (`capacity > count`), so scanning its
  full length would tally reserved, inactive slots (which default to
  node 0 / state S). Pass the active count explicitly.

- counts:

  Optional numeric
  [Column](https://clorton.github.io/razer/reference/Column.md) to
  receive the totals; when omitted an integer result vector is allocated
  and returned instead.

- slot:

  Which slice of `counts` to write when `counts` is supplied; a
  non-negative integer less than `counts`'s slice count. Defaults to
  `0`.

## Value

When `counts` is `NULL`, an integer vector of per-group counts;
otherwise `NULL` invisibly (the result is written into `counts`).

## Details

Two output modes: leave `counts` `NULL` (the default) for an ad-hoc
query and an integer vector of length `n_groups` is allocated and
returned; or pass a numeric
[Column](https://clorton.github.io/razer/reference/Column.md) `counts`
(e.g. a `n_ticks x n_nodes` report) and the totals are written into its
slice `slot` (and `NULL` returned invisibly), avoiding an allocation in
a per-tick model loop.

## See also

[`bincount()`](https://clorton.github.io/razer/reference/bincount.md),
[`bincount_wt()`](https://clorton.github.io/razer/reference/bincount_wt.md),
[`bincount_where_wt()`](https://clorton.github.io/razer/reference/bincount_where_wt.md).

## Examples

``` r
states <- laser_states()
state  <- allocate_scalar("u8",  6L); state$set(c(0, 1, 1, 2, 1, 0))   # S E E I E S
#> NULL
nodeid <- allocate_scalar("u16", 6L); nodeid$set(c(0, 0, 1, 1, 1, 0))
#> NULL
bincount_where(nodeid, 2L, state, "eq", states[["E"]], count = 6L)   # exposed per node: 1 2
#> [1] 1 2
```
