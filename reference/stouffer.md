# Stouffer's intervening-opportunities migration model (Stouffer, 1940).

For each source `i` and destination `j`,
`network[i, j] = k * p_i^a * (p_j / s)^b`, where `s` is the cumulative
population as close or closer to `i` than `j` (excluding the home
population when `include_home = FALSE`); the nearest node (the source)
gets weight 0. Port of laser-core's `stouffer`. The diagonal is 0.

## Usage

``` r
stouffer(pops, distances, k, a, b, include_home)
```

## Arguments

- pops:

  Numeric vector of node populations (length N, non-negative).

- distances:

  Symmetric `N × N` numeric distance matrix.

- k:

  Scaling constant for the flow magnitude (non-negative).

- a:

  Exponent on the origin population.

- b:

  Exponent on the destination/cumulative-population ratio.

- include_home:

  Logical; whether the home population is included in the cumulative
  sum.

## Value

An `N × N` numeric migration-weight matrix (generally not symmetric).
