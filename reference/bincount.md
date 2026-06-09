# Count occurrences of each value, NumPy `bincount`-style, into a buffer.

For each bin `b` in `0..nbins`, counts how many elements of `values`
equal `b` and writes the result into **slice `slot`** of `counts` — the
`nbins` entries `slot*slice_len .. slot*slice_len + nbins` — overwriting
them and leaving the rest untouched. For a scalar `counts` (shape
`(1, n)`) the only slice is `slot = 0`, the whole vector; for a 2-D
report (e.g. `n_ticks x n_nodes` from
[`allocate_vector()`](https://clorton.github.io/razer/reference/allocate_vector.md))
`slot` selects a tick's row. The tally is parallelized with private
per-thread histograms, so there are no write collisions.

## Usage

``` r
bincount(values, nbins, counts, slot = 0L)
```

## Arguments

- values:

  An integer-typed
  [Column](https://clorton.github.io/razer/reference/Column.md) of bin
  indices (`i8`..`u32`), each in `0..nbins`.

- nbins:

  Number of bins; a non-negative integer no greater than `counts`'s
  slice length.

- counts:

  A numeric
  [Column](https://clorton.github.io/razer/reference/Column.md) that
  receives the counts; modified in place.

- slot:

  Which slice of `counts` to write; a non-negative integer less than
  `counts`'s slice count. Defaults to `0`.

## Value

`NULL`, invisibly; the result is written into `counts`.

## See also

[`bincount_wt()`](https://clorton.github.io/razer/reference/bincount_wt.md),
[`bincount_where()`](https://clorton.github.io/razer/reference/bincount_where.md),
[`bincount_where_wt()`](https://clorton.github.io/razer/reference/bincount_where_wt.md).

## Examples

``` r
values <- allocate_scalar("u16", 6L)
values$set(c(0, 1, 1, 2, 2, 2))
#> NULL
counts <- allocate_scalar("i32", 3L)
bincount(values, 3L, counts)
#> NULL
counts$values()   # 1 2 3
#> [1] 1 2 3
```
