# Create a continuous uniform distribution on the half-open interval \[low, high).

Every value in `[low, high)` is equally likely. The mean is
`(low + high) / 2`.

## Usage

``` r
dist_uniform(low, high)
```

## Arguments

- low:

  Inclusive lower bound.

- high:

  Exclusive upper bound; must be strictly greater than `low`.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_uniform(3, 8)   # values in [3, 8), mean 5.5
d$sample_one()
#> [1] 4.855064
```
