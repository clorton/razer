# Fixed-capacity struct-of-arrays population or patch data store.

Mirrors `laser.core.LaserFrame` from Python. Each property occupies a
pre-allocated backing array; `count` tracks how many entries are active.

## Usage

``` r
LaserFrame
```

## Format

An object of class `environment` of length 17.

## Details

**Scalar properties** store one value per entry (length = `capacity`).

**Vector properties** store `length` values per entry, laid out
column-major so that `get_col(name, col)` returns a contiguous memory
slice. This is R's native matrix storage order and makes per-tick
node-state updates fast.

## Methods

### Method `new`

Create a new `LaserFrame`.

#### Arguments

- `capacity`:

  Maximum number of entries (agents or nodes). Must be positive.

- `initial_count`:

  Number of entries active at construction. Pass `-1` (the default) to
  set active count equal to `capacity`.

#### return

A new `LaserFrame` object.

#### export

### Method `count`

Number of currently active entries.

#### export

### Method `capacity`

Total capacity (fixed at construction).

#### export

### Method `scalar_names`

Names of all scalar properties, sorted alphabetically.

#### export

### Method `vector_names`

Names of all vector properties, sorted alphabetically.

#### export

### Method `vector_ncols`

Number of columns in a named vector property.

#### Arguments

- `name`:

  Vector property name.

#### export

### Method `describe`

Human-readable summary of the frame's properties and memory layout.

#### export

### Method `add_scalar_property`

Add a scalar property (one value per entry).

#### Arguments

- `name`:

  Property name. Must not already exist.

- `dtype`:

  `"integer"`, `"real"`, or `"logical"`.

- `default`:

  Fill value for the backing array.

#### export

### Method `add_vector_property`

Add a vector property (`capacity × length`, column-major).

The backing array has `capacity * length` elements. Element
`[entry, col]` is stored at offset `entry + col * capacity`. Column
`col` (all active entries for one time-step) is contiguous.

#### Arguments

- `name`:

  Property name. Must not already exist.

- `length`:

  Number of columns (e.g., `nticks + 1` for a time-series).

- `dtype`:

  `"integer"`, `"real"`, or `"logical"`.

- `default`:

  Fill value.

#### export

### Method `add`

Activate `n` additional entries and return their 1-based index range.

Returns `c(start, end)` (both inclusive, 1-based) so that
`frame$get("prop")[start:end]` addresses the newly activated entries.

#### Arguments

- `n`:

  Number of entries to activate.

#### export

### Method `squash`

Compact active entries, keeping only those where `mask` is `TRUE`.

All scalar and vector properties are squashed in place. `count` is
updated to the number of kept entries.

#### Arguments

- `mask`:

  Logical vector of length `count`.

#### export

### Method `sort_by`

Reorder active scalar properties by `perm` (1-based permutation of
length `count`).

Each property is permuted in parallel across available CPU threads via
Rayon. Only scalar properties are reordered; vector properties
(time-series) are left unchanged.

#### Arguments

- `perm`:

  Integer vector of length `count`. Must be a valid permutation of
  `1:count`.

#### export

### Method `get`

Return the active slice of a scalar property as an R vector.

#### Arguments

- `name`:

  Scalar property name.

#### export

### Method `set`

Overwrite the active slice of a scalar property from an R vector.

Integer and real vectors are accepted for integer and real properties
respectively; integer values are accepted for real properties (coerced).

#### Arguments

- `name`:

  Scalar property name.

- `values`:

  R vector of length `count`.

#### export

### Method `get_col`

Return one column of a vector property as an R vector.

The column is returned for active entries only (length = `count`).
Columns are **1-based**: column 1 is the first column.

#### Arguments

- `name`:

  Vector property name.

- `col`:

  Column index (1-based).

#### export

### Method `set_col`

Overwrite one column of a vector property from an R vector.

#### Arguments

- `name`:

  Vector property name.

- `col`:

  Column index (1-based).

- `values`:

  R vector of length `count`.

#### export

### Method `get_matrix`

Return a vector property as an R matrix of shape `(count, ncols)`.

Because both R and the backing store use column-major layout, each
column of the returned matrix is a direct contiguous copy.

#### Arguments

- `name`:

  Vector property name.

#### export
