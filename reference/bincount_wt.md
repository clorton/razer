# Weighted bincount: sum each element's weight into its bin.

Like
[`bincount()`](https://clorton.github.io/razer/reference/bincount.md),
but accumulates `weights[i]` (rather than 1) into the bin `values[i]`, à
la `numpy.bincount(values, weights = ...)`. The per-bin sums are written
into **slice `slot`** of `counts`. `weights` must be the same length as
`values` and may be any numeric
[Column](https://clorton.github.io/razer/reference/Column.md) (signed,
unsigned, or floating point).

## Usage

``` r
bincount_wt(values, weights, nbins, counts, slot = 0L)
```

## Arguments

- values:

  An integer-typed
  [Column](https://clorton.github.io/razer/reference/Column.md) of bin
  indices (`i8`..`u32`).

- weights:

  A numeric
  [Column](https://clorton.github.io/razer/reference/Column.md), the
  same length as `values`.

- nbins:

  Number of bins; a non-negative integer no greater than `counts`'s
  slice length.

- counts:

  A numeric
  [Column](https://clorton.github.io/razer/reference/Column.md) that
  receives the weighted sums; modified in place.

- slot:

  Which slice of `counts` to write; a non-negative integer less than
  `counts`'s slice count. Defaults to `0`.

## Value

`NULL`, invisibly; the result is written into `counts`.

## See also

[`bincount()`](https://clorton.github.io/razer/reference/bincount.md),
[`bincount_where()`](https://clorton.github.io/razer/reference/bincount_where.md),
[`bincount_where_wt()`](https://clorton.github.io/razer/reference/bincount_where_wt.md).

## Examples

``` r
values  <- allocate_scalar("u16", 5L); values$set(c(0, 0, 1, 2, 2))
#> NULL
weights <- allocate_scalar("f64", 5L); weights$set(c(1.5, 2.5, 4, 1, 3))
#> NULL
counts  <- allocate_scalar("f64", 3L)
bincount_wt(values, weights, 3L, counts)
#> NULL
counts$values()   # 4 4 4
#> [1] 4 4 4
```
