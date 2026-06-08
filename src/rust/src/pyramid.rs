// ════════════════════════════════════════════════════════════════════════════
// pyramid.rs — sampling a discrete distribution by the Vose alias method.
//
// `AliasedDistribution` turns a vector of non-negative integer COUNTS (e.g. the
// population in each age bin of a demographic pyramid) into an O(1)-per-draw sampler
// over the bin indices, weighted by those counts. It is the engine for "build a
// realistic age structure": sample a bin, then (in R) pick a uniform age inside that
// bin's year range.
//
// The Vose alias method (Michael Vose, 1991) precomputes, for each bin, a single
// "alias" partner and an integer threshold so that one draw is just: pick a bin `i`
// uniformly, pick an integer `d` uniformly in [0, total); if `d < probs[i]` keep `i`,
// else jump to `alias[i]`. All arithmetic is INTEGER (no floating-point round-off), a
// faithful port of laser.core's `AliasedDistribution`.
//
// Orientation for readers coming from C / C++ / C# / Python:
//   * `#[extendr]` on the struct generates an opaque external-pointer handle for R;
//     on the `impl` block it exposes each method as `obj$method(...)`. A free `fn`
//     returning the struct is the constructor R calls (`aliased_distribution(...)`).
//   * `Vec<T>` is a growable array (std::vector / list). `&[T]` is a borrowed slice.
//   * `while let (Some(a), Some(b)) = (x.pop(), y.pop())` pops both stacks and runs
//     the body only while BOTH yield a value — `pop()` returns `Option` (nullable).
//   * `par_chunks_mut` (Rayon) splits a mutable slice into contiguous chunks handed
//     to worker threads — the same across-agents parallelism used elsewhere (see the
//     project memory). Each chunk's RNG comes from `rng.rs` (seedable via `set_seed`).
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::rng;

/// A discrete distribution over bin indices `0..n`, sampled by the Vose alias method.
///
/// Construct with [aliased_distribution()] from a vector of non-negative counts; the
/// probability of drawing bin `i` is `counts[i] / sum(counts)`. Draws are 0-based bin
/// indices. The handle is opaque to R.
///
/// @export
#[extendr]
pub struct AliasedDistribution {
    // alias[i] = the partner bin that draws "overflowing" from column i land in
    // (-1 if bin i never overflows — its own mass fills its column exactly).
    alias: Vec<i32>,
    // probs[i] = integer height of bin i's OWN mass within its column, out of `total`.
    // Held as i64 because the construction scales each count by the bin count, which
    // can exceed 32 bits for national-scale populations.
    probs: Vec<i64>,
    // Column height = sum of the original counts; the upper bound of the `d` draw.
    total: i64,
}

impl AliasedDistribution {
    // One draw with a caller-supplied RNG (shared by the single- and batch- paths).
    // `gen_range(0..hi)` is a half-open uniform integer draw, like NumPy's randint.
    fn sample_index<R: Rng + ?Sized>(&self, rng: &mut R) -> i32 {
        let i = rng.gen_range(0..self.alias.len());
        let d = rng.gen_range(0..self.total);
        if d < self.probs[i] { i as i32 } else { self.alias[i] }
    }
}

#[extendr]
impl AliasedDistribution {
    /// Draw a single bin index (0-based) using a thread-local RNG.
    ///
    /// @return A single integer bin index in `0..n_bins()`.
    fn sample_one(&self) -> i32 {
        self.sample_index(&mut rng::single_rng())
    }

    /// Draw `n` bin indices (0-based), returned as an integer vector.
    ///
    /// The draws are split across CPU cores (each with its own thread-local RNG), so
    /// generating a national-scale population in one call is cheap.
    ///
    /// @param n Number of samples to draw; must be non-negative.
    /// @return An integer vector of length `n` of bin indices in `0..n_bins()`.
    fn sample_n(&self, n: i32) -> Vec<i32> {
        assert!(n >= 0, "n must be non-negative, got {n}");
        let n = n as usize;
        let mut out = vec![0i32; n];
        let base = rng::next_call_base();
        out.par_chunks_mut(rng::RNG_CHUNK).enumerate().for_each(|(ci, c)| {
            let mut r = rng::chunk_rng(base, ci); // per-chunk seeded RNG (reproducible)
            for x in c.iter_mut() {
                *x = self.sample_index(&mut r);
            }
        });
        out
    }

    /// The number of bins.
    /// @return An integer.
    fn n_bins(&self) -> i32 {
        self.alias.len() as i32
    }

    /// The total weight (sum of the original counts).
    /// @return A numeric scalar (the sum may exceed R's integer range).
    fn total(&self) -> f64 {
        self.total as f64
    }

    /// The alias table (0-based partner bin per bin, or -1 for none). For inspection.
    /// @return An integer vector of length `n_bins()`.
    fn alias(&self) -> Vec<i32> {
        self.alias.clone()
    }

    /// The per-bin own-mass thresholds (out of `total()`). For inspection.
    /// @return A numeric vector of length `n_bins()`.
    fn probs(&self) -> Vec<f64> {
        self.probs.iter().map(|&p| p as f64).collect()
    }
}

/// Build an [AliasedDistribution] from a vector of non-negative bin counts.
///
/// The probability of drawing bin `i` (0-based) is `counts[i] / sum(counts)`. Counts
/// are rounded to whole numbers; they must be finite and non-negative and sum to a
/// positive total. A typical use is the per-age-bin population of a demographic
/// pyramid (e.g. males + females in each five-year band).
///
/// @param counts A numeric vector of non-negative per-bin counts (length >= 1).
/// @return An `AliasedDistribution` object.
/// @examples
/// d <- aliased_distribution(c(10, 30, 60))   # bin 2 drawn ~60% of the time
/// table(d$sample_n(10000L))                  # 0-based bin indices
/// @export
#[extendr]
fn aliased_distribution(counts: Vec<f64>) -> AliasedDistribution {
    let n = counts.len();
    assert!(n >= 1, "`counts` must have at least one bin");

    // Round/validate the incoming counts to integers.
    let icounts: Vec<i64> = counts
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            assert!(
                c.is_finite() && c >= 0.0,
                "count[{i}] must be finite and non-negative, got {c}"
            );
            c.round() as i64
        })
        .collect();

    let total: i64 = icounts.iter().sum();
    assert!(total > 0, "`counts` must sum to a positive total");

    // Vose setup: scale each count by the bin count so we can compare a bin's mass
    // against the average using only integers — bin i is "small" if probs[i] < total
    // (i.e. count[i] < average), "large" if probs[i] > total. (Derivation: we want
    // count[i] < sum/n  <=>  count[i]*n < sum.)
    let mut probs: Vec<i64> = icounts.iter().map(|&c| c * n as i64).collect();
    let mut alias = vec![-1i32; n];

    // Worklists of under- and over-full bins.
    let mut small: Vec<usize> = (0..n).filter(|&i| probs[i] < total).collect();
    let mut large: Vec<usize> = (0..n).filter(|&i| probs[i] > total).collect();

    // Pair each under-full bin with an over-full one: fill the small column to exactly
    // `total` by borrowing `total - probs[ismall]` from the large bin, then re-classify
    // the large bin (it may now be small, large, or exactly full). Integer arithmetic
    // conserves the total (sum of probs == total*n throughout), so the two worklists
    // drain together — bins left untouched are exactly full and alias to themselves.
    while let (Some(ismall), Some(ilarge)) = (small.pop(), large.pop()) {
        alias[ismall] = ilarge as i32;
        probs[ilarge] -= total - probs[ismall];
        if probs[ilarge] < total {
            small.push(ilarge);
        } else if probs[ilarge] > total {
            large.push(ilarge);
        }
    }

    AliasedDistribution { alias, probs, total }
}

extendr_module! {
    mod pyramid;
    impl AliasedDistribution;
    fn aliased_distribution;
}
