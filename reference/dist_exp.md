# Create an exponential distribution with rate `rate`.

Draws are strictly positive with mean `1 / rate` and variance
`1 / rate^2`.

## Usage

``` r
dist_exp(rate)
```

## Arguments

- rate:

  Rate parameter λ; must be finite and positive.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_exp(0.5)   # mean 2
d$sample_one()
#> [1] 3.925157
```
