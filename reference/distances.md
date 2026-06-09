# Great-circle distance matrix between geographic points.

Ports laser-core's `distance` (all-pairs case): given the latitudes and
longitudes of `N` points in decimal degrees, returns the symmetric
`N × N` matrix whose `[i, j]` entry is the haversine great-circle
distance, in **kilometres**, between point `i` and point `j`. The
diagonal is zero.

## Usage

``` r
distances(latitude, longitude)
```

## Arguments

- latitude:

  Numeric vector of latitudes in decimal degrees, in `[-90, 90]`.

- longitude:

  Numeric vector of longitudes in decimal degrees, in `[-180, 180]`.
  Must be the same length as `latitude`.

## Value

An `N × N` numeric matrix of pairwise distances in kilometres
(symmetric, zero diagonal), where `N` is the number of points.

## Details

The haversine formula (see
<https://en.wikipedia.org/wiki/Haversine_formula>) is evaluated with a
mean Earth radius of 6371 km, matching laser-core.

## Examples

``` r
# London and Paris, ~344 km apart:
d <- distances(c(51.5074, 48.8566), c(-0.1278, 2.3522))
d[1, 2]
#> [1] 343.5561
```
