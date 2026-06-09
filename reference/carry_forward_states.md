# Carry census counters forward, and optionally total some of them.

For each Column in `carry`, calls
[`carry_forward()`](https://clorton.github.io/razer/reference/carry_forward.md)
(copying column `tick` onto `tick + 1`). If `total` is supplied, then
after carrying it sets `total`'s column `tick + 1` to the elementwise
sum of the `summands` Columns at `tick + 1` — for example carrying `S`,
`I`, `R` forward and totalling them into `N` so the current per-node
population is available to
[`calc_foi()`](https://clorton.github.io/razer/reference/calc_foi.md)
(and stays correct as births, deaths, and imports change the states).

## Usage

``` r
carry_forward_states(carry, tick, total = NULL, summands = carry)
```

## Arguments

- carry:

  A list of 2-D census
  [Column](https://clorton.github.io/razer/reference/Column.md)s to
  carry forward.

- tick:

  0-based source tick; column `tick` is copied onto `tick + 1`.

- total:

  Optional [Column](https://clorton.github.io/razer/reference/Column.md)
  to receive the running total at column `tick + 1`.

- summands:

  List of Columns to sum into `total` (defaults to `carry`).

## Value

`NULL`, invisibly; the Columns are modified in place.

## Examples

``` r
if (FALSE) { # \dontrun{
# carry S, I, R forward and keep N = S + I + R up to date:
carry_forward_states(list(nodes$S, nodes$I, nodes$R), tick, total = nodes$N)
} # }
```
