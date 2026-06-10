# Create a beta distribution on the open interval (0, 1).

Parameterized by the two positive shape parameters α and β. The mean is
`alpha / (alpha + beta)`.

## Usage

``` r
dist_beta(alpha, beta)
```

## Arguments

- alpha:

  First shape parameter (α); must be finite and positive.

- beta:

  Second shape parameter (β); must be finite and positive.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_beta(2, 5)   # mean 2/7 ≈ 0.286, support (0, 1)
d$sample_one()
#> [1] 0.2400854
```
