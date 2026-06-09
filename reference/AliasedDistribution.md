# A discrete distribution over bin indices `0..n`, sampled by the Vose alias method.

Construct with
[`aliased_distribution()`](https://clorton.github.io/razer/reference/aliased_distribution.md)
from a vector of non-negative counts; the probability of drawing bin `i`
is `counts[i] / sum(counts)`. Draws are 0-based bin indices. The handle
is opaque to R.

## Usage

``` r
AliasedDistribution
```

## Format

An object of class `environment` of length 6.

## Methods

### Method `sample_one`

Draw a single bin index (0-based) using a thread-local RNG.

#### return

A single integer bin index in `0..n_bins()`.

### Method `sample_n`

Draw `n` bin indices (0-based), returned as an integer vector.

The draws are split across CPU cores (each with its own thread-local
RNG), so generating a national-scale population in one call is cheap.

#### Arguments

- `n`:

  Number of samples to draw; must be non-negative.

#### return

An integer vector of length `n` of bin indices in `0..n_bins()`.

### Method `n_bins`

The number of bins.

#### return

An integer.

### Method `total`

The total weight (sum of the original counts).

#### return

A numeric scalar (the sum may exceed R's integer range).

### Method `alias`

The alias table (0-based partner bin per bin, or -1 for none). For
inspection.

#### return

An integer vector of length `n_bins()`.

### Method `probs`

The per-bin own-mass thresholds (out of `total()`). For inspection.

#### return

A numeric vector of length `n_bins()`.
