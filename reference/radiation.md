# Radiation migration-network model (Simini et al., Nature 2012).

For each source `i`, destinations are ranked by distance and the weight
to destination `j` is `k * p_i * p_j / (p_i + s) / (p_i + p_j + s)`,
where `s` is the total population of all nodes as close or closer to `i`
than `j` (excluding the home population `p_i` when
`include_home = FALSE`). Port of laser-core's `radiation`. The diagonal
is 0.

## Usage

``` r
radiation(pops, distances, k, include_home)
```

## Arguments

- pops:

  Numeric vector of node populations (length N, non-negative).

- distances:

  Symmetric `N × N` numeric distance matrix.

- k:

  Scaling constant for the flow magnitude (non-negative).

- include_home:

  Logical; whether the home (source) population is included in the "as
  close or closer" cumulative sum.

## Value

An `N × N` numeric migration-weight matrix (generally not symmetric).
