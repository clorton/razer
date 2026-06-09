# Build a values map: broadcast a value into a `n_ticks x n_nodes` grid [Column](https://clorton.github.io/razer/reference/Column.md).

Expands a flexible per-time and/or per-node value into a full
`n_ticks x n_nodes` f64
[Column](https://clorton.github.io/razer/reference/Column.md)
(slice-major: each tick's per-node row contiguous), suitable for passing
to [`calc_foi()`](https://clorton.github.io/razer/reference/calc_foi.md)
as `beta` or `seasonality`. This is razer's equivalent of LASER's
*ValuesMap*. The shape of `value` selects how it is broadcast:

## Usage

``` r
values_map(value, n_ticks, n_nodes)
```

## Arguments

- value:

  A scalar, a length-`n_nodes` or length-`n_ticks` numeric vector, or a
  `n_ticks x n_nodes` numeric matrix.

- n_ticks:

  Number of time slices (rows of the grid).

- n_nodes:

  Number of nodes (columns of the grid).

## Value

An f64 [Column](https://clorton.github.io/razer/reference/Column.md) of
shape `n_ticks x n_nodes`.

## Details

- **scalar** (length 1) — constant over time and space.

- **length `n_nodes`** — varies by node, constant over time (per-node).

- **length `n_ticks`** — varies by time, constant over space (per-tick).

- **`n_ticks x n_nodes` matrix** — varies by both; used as-is.

When `n_ticks == n_nodes` a bare vector is ambiguous; it is treated as
per-node. Pass an explicit matrix to vary by both in that case.

## Examples

``` r
g <- values_map(0.5, 10L, 3L)        # constant 0.5 everywhere
dim(g$values())                       # 10 3
#> [1] 10  3
values_map(c(1, 2, 3), 10L, 3L)       # per-node (length n_nodes)
#> <pointer: 0x56192a773ac0>
#> attr(,"class")
#> [1] "Column"
values_map(seq_len(10L), 10L, 3L)     # per-tick (length n_ticks)
#> <pointer: 0x56192a754d00>
#> attr(,"class")
#> [1] "Column"
```
