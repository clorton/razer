# Gravity migration-network model.

`network[i, j] = k * pops[i]^a * pops[j]^b / distance[i, j]^c`, with a
zero diagonal (no self-migration). Port of laser-core's `gravity`.

## Usage

``` r
gravity(pops, distances, k, a, b, c)
```

## Arguments

- pops:

  Numeric vector of node populations (length N, non-negative).

- distances:

  Symmetric `N × N` numeric distance matrix (e.g. from
  [`distances()`](https://clorton.github.io/razer/reference/distances.md)).

- k:

  Scaling constant for the overall flow magnitude (non-negative).

- a:

  Exponent on the origin population.

- b:

  Exponent on the destination population.

- c:

  Exponent on the distance (larger `c` = stronger distance decay).

## Value

An `N × N` numeric matrix; `[i, j]` is the migration weight from node i
to node j.
