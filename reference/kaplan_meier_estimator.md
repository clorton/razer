# Build a [KaplanMeierEstimator](https://clorton.github.io/razer/reference/KaplanMeierEstimator.md) from cumulative deaths by year.

`cumulative_deaths[y]` is the number of a synthetic cohort dead by the
end of year `y`; it must be non-negative and monotonically
non-decreasing. Values are rounded to whole numbers. (A leading zero is
prepended internally; do not include it yourself.)

## Usage

``` r
kaplan_meier_estimator(cumulative_deaths)
```

## Arguments

- cumulative_deaths:

  A non-decreasing numeric vector of cumulative deaths by year (length
  \>= 1).

## Value

A `KaplanMeierEstimator` object.

## Examples

``` r
# toy life table: 10 deaths/year for 80 years, then 100/year for 21 more
km <- kaplan_meier_estimator(cumsum(c(rep(10, 80), rep(100, 21))))
```
