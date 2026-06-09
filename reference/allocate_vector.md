# Allocate a fresh, zero-filled 2-D property array (a per-slot report buffer).

Returns an opaque
[Column](https://clorton.github.io/razer/reference/Column.md) of
`n_slices * slice_len` elements laid out SLICE-MAJOR (row-major):
`n_slices` slices, each a contiguous run of `slice_len` elements. Slice
`s` is the block `s*slice_len .. (s+1)*slice_len`, so indexing the FIRST
dimension yields a contiguous array. The conventional use is a
time-series report with the **first dimension time and the second node**
— `allocate_vector(dtype, n_ticks, n_nodes)` — so each tick's per-node
values are contiguous (cache-friendly for the step kernels that fill one
tick at a time). `$values()` reads it back as an `n_slices × slice_len`
(e.g. `n_ticks × n_nodes`) R matrix, so row `t` is tick `t`'s vector.

## Usage

``` r
allocate_vector(dtype, n_slices, slice_len)
```

## Arguments

- dtype:

  Element type (see
  [`allocate_scalar()`](https://clorton.github.io/razer/reference/allocate_scalar.md)
  for the accepted names).

- n_slices:

  Number of slices — the first/outer dimension (e.g. the tick count). A
  non-negative integer.

- slice_len:

  Contiguous length of each slice — the inner dimension (e.g. the node
  count). A non-negative integer.

## Value

A [Column](https://clorton.github.io/razer/reference/Column.md) of shape
`n_slices × slice_len`, all elements zero.

## Examples

``` r
recoveries <- allocate_vector("u32", 4L, 3L)   # 4 ticks x 3 nodes
dim(recoveries$values())                        # 4 3
#> [1] 4 3
```
