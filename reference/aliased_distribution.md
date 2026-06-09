# Build an [AliasedDistribution](https://clorton.github.io/razer/reference/AliasedDistribution.md) from a vector of non-negative bin counts.

The probability of drawing bin `i` (0-based) is
`counts[i] / sum(counts)`. Counts are rounded to whole numbers; they
must be finite and non-negative and sum to a positive total. A typical
use is the per-age-bin population of a demographic pyramid (e.g. males +
females in each five-year band).

## Usage

``` r
aliased_distribution(counts)
```

## Arguments

- counts:

  A numeric vector of non-negative per-bin counts (length \>= 1).

## Value

An `AliasedDistribution` object.

## Examples

``` r
d <- aliased_distribution(c(10, 30, 60))   # bin 2 drawn ~60% of the time
table(d$sample_n(10000L))                  # 0-based bin indices
#> 
#>    0    1    2 
#>  946 3025 6029 
```
