# Create a gamma distribution parameterized by shape and scale.

Uses the shape–scale (k, θ) parameterization: the mean is
`shape * scale` and the variance is `shape * scale^2`. Draws are
strictly positive, which makes the gamma a natural choice for
right-skewed, always-positive durations.

## Usage

``` r
dist_gamma(shape, scale)
```

## Arguments

- shape:

  Shape parameter k; must be finite and positive.

- scale:

  Scale parameter θ; must be finite and positive.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_gamma(2, 3)   # mean 6, variance 18
d$sample_one()
#> [1] 6.242076
```
