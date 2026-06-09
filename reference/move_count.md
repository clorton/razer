# Apply a per-node transition count to the census at `tick + 1`.

Subtracts `counts` from the `from` state and adds it to the `to` state
at census column `tick + 1` (the working column the model has already
carried forward). Either side may be `NULL` to skip it — e.g. a death is
a one-sided decrement (`to = NULL`) and a birth a one-sided increment
(`from = NULL`).

## Usage

``` r
move_count(from, to, counts, tick)
```

## Arguments

- from:

  A 2-D census
  [Column](https://clorton.github.io/razer/reference/Column.md) to
  decrement, or `NULL`.

- to:

  A 2-D census
  [Column](https://clorton.github.io/razer/reference/Column.md) to
  increment, or `NULL`.

- counts:

  Integer vector of per-node counts (length `n_nodes`).

- tick:

  0-based source tick; the delta is applied at column `tick + 1`.

## Value

`NULL`, invisibly; the Columns are modified in place.

## Examples

``` r
if (FALSE) { # \dontrun{
inf <- transmission(state, timer, nodeid, count, nodes$foi, t, states[["E"]], inc_dur)
move_count(nodes$S, nodes$E, inf, t)   # S -> E
} # }
```
