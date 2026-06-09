# Create a normal (Gaussian) distribution.

The second argument is the **variance** (σ²), not the standard
deviation, to match the way variance is usually quoted in statistical
models. The standard deviation passed to the underlying sampler is
`sqrt(variance)`.

## Usage

``` r
dist_normal(mean, variance)
```

## Arguments

- mean:

  Mean (μ) of the distribution.

- variance:

  Variance (σ²); must be non-negative.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_normal(7, 4)   # mean 7, variance 4 (sd 2)
d$sample_one()
#> [1] 7.290824
```
