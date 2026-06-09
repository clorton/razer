# Load a population-pyramid CSV into a numeric matrix.

The file is expected to have the schema used by laser.core pyramids:

## Usage

``` r
load_pyramid_csv(file)
```

## Arguments

- file:

  Path to the CSV file.

## Value

An integer matrix with columns `start`, `end`, `M`, `F` (one row per
band).

## Details

- a header line exactly `"Age,M,F"`;

- then one line per age band `"low-high,males,females"` (all
  non-negative integers);

- a final open-ended band `"max+,males,females"`.

The returned matrix has one row per band and the integer columns
`start`, `end`, `M`, `F`. The final `"max+"` band is stored as a
single-year bucket (`end == start`), matching laser.core.

## Errors

Stops if the header is not `"Age,M,F"`, if any data line is malformed,
or if the start/end ages are not strictly ascending.

## Examples

``` r
if (FALSE) { # \dontrun{
pyramid <- load_pyramid_csv("USA_pyramid_2020.csv")
} # }
```
