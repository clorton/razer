# Create a degenerate (constant) distribution that always returns `value`.

Use this as a fixed-duration drop-in wherever a `Distribution` is
required — e.g. a deterministic infectious period of exactly `value`
ticks.

## Usage

``` r
dist_constant(value)
```

## Arguments

- value:

  The constant value returned by every draw.

## Value

A `Distribution` object.

## Examples

``` r
d <- dist_constant(10)   # always 10
d$sample_one()
#> [1] 10
```
