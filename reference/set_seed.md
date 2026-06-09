# Set the global random seed, making subsequent razer runs reproducible.

After `set_seed(s)`, every kernel's randomness is a deterministic
function of `s` and the order of kernel calls — identical on every run
and every machine, regardless of CPU/thread count. Call it once at the
start of a script.
[`unset_seed()`](https://clorton.github.io/razer/reference/unset_seed.md)
reverts to a fresh (entropy-seeded) RNG.

## Usage

``` r
set_seed(seed)
```

## Arguments

- seed:

  A finite, non-negative number (used as a 64-bit seed).

## Value

`NULL`, invisibly.

## Examples

``` r
set_seed(42)
#> NULL
```
