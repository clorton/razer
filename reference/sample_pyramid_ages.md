# Sample realistic per-agent ages (in days) from a population pyramid.

Builds an
[`aliased_distribution()`](https://clorton.github.io/razer/reference/aliased_distribution.md)
over the per-band populations (`M + F`), draws a band for each agent in
proportion to its population, then a uniform day within that band's year
range `[start, end + 1)` — so an agent in the `0-4` band gets an age
uniformly in `[0, 5)` years.

## Usage

``` r
sample_pyramid_ages(pyramid, n)
```

## Arguments

- pyramid:

  An integer matrix with columns `start`, `end`, `M`, `F`, as returned
  by
  [`load_pyramid_csv()`](https://clorton.github.io/razer/reference/load_pyramid_csv.md).

- n:

  Number of agent ages to draw.

## Value

An integer vector of length `n` of ages in whole days.

## Details

Band selection uses the package's internal (thread-local, not
R-seedable) RNG via the alias sampler; the within-band day uses R's RNG
(`set.seed`-able).

## Examples

``` r
pyramid <- rbind(c(0, 4, 100, 98), c(5, 9, 90, 92), c(10, 10, 5, 6))
colnames(pyramid) <- c("start", "end", "M", "F")
ages_days <- sample_pyramid_ages(pyramid, 1000L)
```
