# Allocate a fresh, zero-filled property array of a given type and length.

Returns an opaque
[Column](https://clorton.github.io/razer/reference/Column.md) handle
backed by a Rust-owned buffer. Choose the narrowest type that fits the
data to minimize memory (e.g. `"u8"` for a small set of disease-state
codes).

## Usage

``` r
allocate_scalar(dtype, count)
```

## Arguments

- dtype:

  Element type, one of `"i8"`, `"u8"`, `"i16"`, `"u16"`, `"i32"`,
  `"u32"`, `"f32"`, `"f64"` (aliases: `"int8"`, `"uint8"`/`"raw"`, …,
  `"integer"` = i32, `"double"`/`"real"` = f64, `"single"` = f32).

- count:

  Array length (number of elements); a non-negative integer.

## Value

A [Column](https://clorton.github.io/razer/reference/Column.md) object
whose elements are all zero.

## Examples

``` r
state <- allocate_scalar("u8", 5L)
state$dtype()    # "u8"
#> [1] "u8"
state$length()   # 5
#> [1] 5
state$values()   # 0 0 0 0 0
#> [1] 0 0 0 0 0
```
