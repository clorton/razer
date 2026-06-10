# A Kaplan–Meier sampler over a cumulative-deaths-by-year life table.

Construct with
[`kaplan_meier_estimator()`](https://clorton.github.io/razer/reference/kaplan_meier_estimator.md)
from a non-decreasing vector of cumulative deaths by year. The handle is
opaque to R.

## Usage

``` r
KaplanMeierEstimator
```

## Format

An object of class `environment` of length 3.

## Methods

### Method `predict_year_of_death`

Predict a year of death for each individual given their current age in
years.

For each age, samples a year of death `>= age` and `<= max_year`,
conditioned on survival to that age (Kaplan–Meier). Ages must be in
`0..=max_year`.

#### Arguments

- `ages_years`:

  Integer vector of current ages in whole years.

- `max_year`:

  Maximum year of death to consider; pass a negative value (e.g. `-1L`)
  to use the last year in the life table.

#### return

An integer vector of predicted years of death (same length as input).

### Method `predict_age_at_death`

Predict an age at death (in DAYS) for each individual given their age in
days.

Samples the year of death as in `predict_year_of_death()`, then a day
within that year: a uniform day of a later year, or — if death falls in
the individual's current year — a uniform day at or after their current
day-of-year (so the predicted age at death is never earlier than the
current age). Ages in days must be `< (max_year + 1) * 365`.

#### Arguments

- `ages_days`:

  Integer vector of current ages in whole days.

- `max_year`:

  Maximum year of death to consider; pass a negative value (e.g. `-1L`)
  to use the last year in the life table.

#### return

An integer vector of predicted ages at death in days (same length).

### Method `cumulative_deaths`

The cumulative-deaths-by-year table (without the internal leading zero).

#### return

A numeric vector of length equal to the number of years.

## Examples

``` r
## ---- Method `predict_year_of_death` ---- ##
km <- kaplan_meier_estimator(cumsum(c(rep(10, 80), rep(100, 21))))
km$predict_year_of_death(c(40L, 50L, 60L), -1L)
#> [1] 82 82 83

## ---- Method `predict_age_at_death` ---- ##
km <- kaplan_meier_estimator(cumsum(c(rep(10, 80), rep(100, 21))))
km$predict_age_at_death(c(40L, 50L, 60L) * 365L, -1L)
#> [1] 36369 31684 31161

```
