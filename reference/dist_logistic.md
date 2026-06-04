# Create a logistic distribution with the given location and scale.

Symmetric about `location` (its mean and median); the variance is
`scale^2 * pi^2 / 3`.

## Usage

``` r
dist_logistic(location, scale)
```

## Arguments

- location:

  Location parameter μ (the mean); must be finite.

- scale:

  Scale parameter s; must be finite and positive.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_logistic(4, 2)   # mean 4
d$sample_one()
#> [1] 3.316977
```
