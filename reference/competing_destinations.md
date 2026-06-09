# Competing-destinations migration model (Fotheringham, 1984).

Starts from the
[`gravity()`](https://clorton.github.io/razer/reference/gravity.md)
weights and multiplies each `[i, j]` by an accessibility term
`(Σ_k p_k^b / d_jk^c)^delta`, summed over competing destinations
`k ∉ {i, j}`. Port of laser-core's `competing_destinations`. The
diagonal is 0.

## Usage

``` r
competing_destinations(pops, distances, k, a, b, c, delta)
```

## Arguments

- pops:

  Numeric vector of node populations (length N, non-negative).

- distances:

  Symmetric `N × N` numeric distance matrix.

- k:

  Scaling constant for the flow magnitude (non-negative).

- a:

  Exponent on the origin population (gravity term).

- b:

  Exponent on the destination population (gravity and competition
  terms).

- c:

  Exponent on distance (gravity and competition terms).

- delta:

  Exponent on the competing-destinations accessibility term.

## Value

An `N × N` numeric migration-weight matrix (generally not symmetric).
