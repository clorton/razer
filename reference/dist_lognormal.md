# Create a log-normal distribution.

A variable whose natural logarithm is `Normal(meanlog, sdlog)`. Draws
are strictly positive. The median is `exp(meanlog)` and the mean is
`exp(meanlog + sdlog^2 / 2)`. `meanlog` and `sdlog` are the log-space
parameters, matching R's `qlnorm(p, meanlog, sdlog)`.

## Usage

``` r
dist_lognormal(meanlog, sdlog)
```

## Arguments

- meanlog:

  Mean of the underlying normal (in log space); must be finite.

- sdlog:

  Standard deviation of the underlying normal; must be finite and
  non-negative.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_lognormal(0, 0.5)   # median 1
d$sample_one()
#> [1] 1.402959
```
