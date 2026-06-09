# A Rust-owned, dtype-tagged property array (1-D scalar, or a 2-D vector report).

Allocate one with
[`allocate_scalar()`](https://clorton.github.io/razer/reference/allocate_scalar.md)
(`nrows` elements, `ncols == 1`) or
[`allocate_vector()`](https://clorton.github.io/razer/reference/allocate_vector.md)
(`nrows`-per-slot × `ncols` slots, stored COLUMN-MAJOR so a whole
slot/time-slice is contiguous). The data is held in Rust and exposed to
R only as an opaque handle; use `$values()` to copy a snapshot back into
an R vector (or matrix, when `ncols > 1`) for inspection, `$fill()` /
`$set()` to write, and `$length()` / `$dtype()` to query. The simulation
step kernels operate on the buffer in place with no copies.

## Usage

``` r
Column
```

## Format

An object of class `environment` of length 8.

## Methods

### Method `length`

Number of elements in the array.

#### return

An integer length.

### Method `dtype`

The element data type as a string (e.g. `"u8"`, `"f32"`).

#### return

A length-1 character vector.

### Method `values`

Copy the array into an R vector for inspection (NOT a view — a
snapshot).

Integer-width types (i8, u8, i16, u16, i32) widen to R `integer`; `u32`,
`f32`, and `f64` widen to R `double` (since `u32` overflows R's signed
32-bit integer). This O(n) copy is the only place data leaves Rust.

For a 2-D column (`n_slices > 1`, from
[`allocate_vector()`](https://clorton.github.io/razer/reference/allocate_vector.md))
the result carries a `dim` attribute and reads back as an
`n_slices × slice_len` R matrix (e.g. `n_ticks × n_nodes`), so row `t`
is tick `t`'s per-node vector; otherwise a plain vector. The snapshot is
transposed during the copy (our buffer is slice-major, R matrices are
column-major) — inexpensive, inspection-only.

#### return

A numeric vector (or matrix) — integer or double — of
[`length()`](https://rdrr.io/r/base/length.html) elements.

### Method `fill`

Set every element to `value`, cast to the array's data type.

For integer-typed arrays the value is truncated toward zero (e.g. `2.9`
becomes `2`); out-of-range values wrap per Rust's `as` cast.

#### Arguments

- `value`:

  A single numeric value to broadcast across the array.

### Method `set`

Overwrite the array from an R numeric vector (integer or double).

The input length must equal
[`length()`](https://rdrr.io/r/base/length.html). Each element is cast
to the array's data type (integer-typed arrays truncate toward zero).
Useful for setup — e.g. writing per-agent node assignments or seeding
initial states from R.

#### Arguments

- `values`:

  A numeric vector of length
  [`length()`](https://rdrr.io/r/base/length.html).

### Method `col`

Read one column (`slot`) of a 2-D Column as an R vector snapshot.

Returns the `slice_len` values in column `slot` (e.g. all nodes for one
tick), widened to R `integer`/`double` like `values()`. For a scalar
column the only valid `slot` is 0 (the whole vector).

#### Arguments

- `slot`:

  0-based column index, less than the number of columns.

#### return

A numeric vector of length `slice_len`.

### Method `set_col`

Write `values` into one column (`slot`) of a 2-D Column, in place.

`values` must have length `slice_len`; each element is cast to the
column's data type. Lets the caller update a single column (e.g. one
tick's per-node slice — a derived population total) without rewriting
the whole buffer.

#### Arguments

- `slot`:

  0-based column index, less than the number of columns.

- `values`:

  A numeric vector of length `slice_len`.

### Method `squash`

Compact the first `length(keep)` elements in place, keeping those
flagged `TRUE`.

Drops elements where `keep` is `FALSE`/`NA`, shifting the survivors to
the front (order preserved), and returns the number kept. Use it to
reclaim the slots of deceased agents: apply the SAME `keep` mask to
every per-agent Column (so they stay aligned) and set the active count
to the returned value. The R
[`squash()`](https://clorton.github.io/razer/reference/squash.md) helper
does exactly this across a people environment. Only valid for a 1-D
(scalar) Column.

#### Arguments

- `keep`:

  A logical vector whose length is at most the column length (typically
  the active agent count).

#### return

The number of kept elements (an integer); elements past it are left
as-is.
