# A parameterized probability distribution that can be sampled repeatedly.

Build one with a family constructor such as
[`dist_normal()`](https://clorton.github.io/razer/reference/dist_normal.md)
or
[`dist_constant()`](https://clorton.github.io/razer/reference/dist_constant.md),
then either sample it from Rust via
[`sample()`](https://rdrr.io/r/base/sample.html) (Pattern B: the caller
owns the RNG and passes it in), or from R via `$sample_one()` /
`$sample_n()`.

## Usage

``` r
Distribution
```

## Format

An object of class `environment` of length 2.

## Details

Draws are always floating-point. Callers that need integer values (e.g.
a whole-tick state timer) are responsible for rounding or truncating as
appropriate — the simulation kernels round to the nearest tick.

The handle is opaque to R — it is passed to simulation kernels (e.g.
`step_sir` or `transmission`) by reference, so the same object can be
reused every tick and shared across all worker threads.

## Methods

### Method `sample_one`

Draw a single sample using a thread-local RNG.

Convenience for interactive use. Simulation kernels do not call this —
they use the internal sampler with an explicit, reusable RNG for
performance.

#### return

A single numeric (double) draw from the distribution.

### Method `sample_n`

Draw `n` samples using a thread-local RNG, returned as a numeric vector.

Drawing a whole batch in one call avoids per-sample R↔Rust overhead,
which makes it practical to validate the sampler against large empirical
samples (e.g. one million draws) from R.

#### Arguments

- `n`:

  Number of samples to draw; must be non-negative.

#### return

A numeric (double) vector of length `n`.
