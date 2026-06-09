# Carry a per-node counter forward one tick: copy column `tick` onto `tick + 1`.

For a 2-D report
[Column](https://clorton.github.io/razer/reference/Column.md) (e.g.
`n_ticks+1 × n_nodes`), copies the contiguous slice for `tick` onto the
slice for `tick + 1`. This seeds the next tick with the current counts
so a dynamics kernel can then update it in place — keeping the census
invariant `count[t+1] = count[t] ± delta`. Call it once per state that
must persist across ticks (e.g. S, I, R for an SIR model; add E for
SEIR, or a user-defined "V" vaccinated count). Works for any element
type.

## Usage

``` r
carry_forward(counter, tick)
```

## Arguments

- counter:

  A 2-D Column to carry forward.

- tick:

  0-based source tick; column `tick` is copied onto column `tick+1`.

## Value

`NULL` (invisibly); `counter` is modified in place.

## Examples

``` r
counts <- allocate_vector("i32", 3L, 2L)   # 3 ticks x 2 nodes
counts$set(c(5L, 7L, 0L, 0L, 0L, 0L))      # tick 0 = (5, 7)
#> NULL
carry_forward(counts, 0L)
#> NULL
counts$values()[2L, ]                       # 5 7  (carried to tick 1)
#> [1] 5 7
```
