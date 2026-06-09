# Cap each row sum of a network matrix at `max_rowsum`.

Rows whose total exceeds `max_rowsum` are scaled down proportionally so
they sum to exactly `max_rowsum`; smaller rows are left unchanged. Port
of laser-core's `row_normalizer`. Useful to bound the fraction of a
node's force of infection that is exported before passing the matrix as
the transmission `network`.

## Usage

``` r
row_normalizer(network, max_rowsum)
```

## Arguments

- network:

  A square non-negative numeric matrix.

- max_rowsum:

  Maximum allowed row sum, in `[0, 1]`.

## Value

The row-capped matrix, same shape as `network`.
