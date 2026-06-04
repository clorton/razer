# Create a Poisson distribution with rate (mean) `lambda`.

Draws are non-negative integer counts (returned as doubles) with mean
and variance both equal to `lambda`. Useful for count-valued durations.

## Usage

``` r
dist_poisson(lambda)
```

## Arguments

- lambda:

  Rate / mean λ; must be finite and positive.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_poisson(5)   # mean 5, integer-valued draws
d$sample_one()
#> [1] 6
```
