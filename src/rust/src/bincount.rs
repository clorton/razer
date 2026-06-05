// ════════════════════════════════════════════════════════════════════════════
// bincount.rs — a parallel, collision-free histogram over a Column, à la NumPy's
// numpy.bincount().
//
// `bincount(values, nbins, counts)` computes, for each bin b in 0..nbins, how many
// elements of `values` equal b, and writes the result into the caller-provided
// `counts` Column. `values` is any INTEGER-typed Column (its elements are bin
// indices); `counts` is any numeric Column with at least `nbins` elements.
//
// Parallelism strategy (the reason this isn't just `for v in values { counts[v]+=1 }`):
// a shared `counts[v] += 1` from many threads would race on the same bin. Instead
// each Rayon task accumulates into its OWN private histogram buffer (no sharing,
// no locks, no collisions); we then ZERO the used range of `counts` and fold every
// task's buffer into it. Zeroing first is required because that final step
// ACCUMULATES (`counts[b] += local[b]`) rather than overwriting.
//
// Orientation for readers coming from C / C++ / C# / Python:
//   * Generics + trait bounds (`fn histogram<V: BinIndex>`) are like C++ templates
//     with concepts: ONE implementation is compiled (monomorphized) for each
//     element type, so we don't hand-write a copy per int width. The `match` on the
//     Column's variant is only to recover the concrete element type from the tagged
//     union before calling the single generic kernel.
//   * `trait BinIndex { fn to_index(self) -> usize; }` is an interface; the
//     `impl_bin_index!` macro_rules below stamps out the trivial impl for each
//     numeric type at compile time (a hygienic code generator, not a runtime call).
//   * `par_iter().fold(init, f).reduce(init, g)` is Rayon's map-reduce: `fold`
//     builds one accumulator PER TASK (init = a zeroed histogram), `reduce` combines
//     them pairwise. Think OpenMP `reduction` or .NET `Parallel.ForEach` with
//     thread-local state.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use crate::column::{Column, Storage};

// ── values: anything usable as a bin index ──────────────────────────────────────

// `Sync` lets a shared `&[V]` be read from many Rayon threads at once.
trait BinIndex: Copy + Sync {
    fn to_index(self) -> usize;
}

// Stamp out `impl BinIndex` for each integer width. Negative signed values wrap to
// a huge `usize` and are then caught by the bounds check on the per-task buffer
// (an out-of-range index panics, surfacing as an R error) — bin indices must be in
// `0..nbins`.
macro_rules! impl_bin_index {
    ($($t:ty),*) => { $(
        impl BinIndex for $t {
            #[inline]
            fn to_index(self) -> usize { self as usize }
        }
    )* };
}
impl_bin_index!(i8, u8, i16, u16, i32, u32);

// Build one private histogram per Rayon task, then reduce them to a single
// per-bin total vector. No bin is ever written by two threads at once.
fn histogram<V: BinIndex>(values: &[V], nbins: usize) -> Vec<u64> {
    values
        .par_iter()
        .fold(
            || vec![0u64; nbins],                       // per-task zeroed buffer
            |mut local, &v| {
                local[v.to_index()] += 1;               // collision-free local bump
                local
            },
        )
        .reduce(
            || vec![0u64; nbins],
            |mut a, b| {                                // combine two task buffers
                for i in 0..nbins {
                    a[i] += b[i];
                }
                a
            },
        )
}

// ── counts: any numeric output buffer ───────────────────────────────────────────

trait BinCount: Copy {
    fn zero() -> Self;
    fn add_u64(&mut self, x: u64);   // unweighted totals
    fn add_f64(&mut self, x: f64);   // weighted totals
}

macro_rules! impl_bin_count {
    ($($t:ty),*) => { $(
        impl BinCount for $t {
            #[inline]
            fn zero() -> Self { 0 as $t }
            #[inline]
            fn add_u64(&mut self, x: u64) { *self += x as $t; }
            #[inline]
            fn add_f64(&mut self, x: f64) { *self += x as $t; }
        }
    )* };
}
impl_bin_count!(i8, u8, i16, u16, i32, u32, f32, f64);

// Zero `counts[0..nbins]`, then accumulate the integer per-bin totals into it.
fn write_counts_u64<C: BinCount>(counts: &mut [C], totals: &[u64], nbins: usize) {
    for b in 0..nbins {
        counts[b] = C::zero();      // required: the next step ADDS into counts
    }
    for b in 0..nbins {
        counts[b].add_u64(totals[b]);
    }
}

// As above, for floating-point (weighted) per-bin totals.
fn write_counts_f64<C: BinCount>(counts: &mut [C], totals: &[f64], nbins: usize) {
    for b in 0..nbins {
        counts[b] = C::zero();
    }
    for b in 0..nbins {
        counts[b].add_f64(totals[b]);
    }
}

/// Count occurrences of each value, NumPy `bincount`-style, into a buffer.
///
/// For each bin `b` in `0..nbins`, computes how many elements of `values` equal
/// `b` and writes that into `counts[b]`. `values` must be an integer-typed
/// [Column] whose elements lie in `0..nbins` (they are used as bin indices);
/// `counts` must be a numeric [Column] with at least `nbins` elements and is
/// **overwritten** in its first `nbins` entries (entries at or beyond `nbins` are
/// left untouched). The work is parallelized with private per-thread histograms,
/// so there are no write collisions.
///
/// @param values An integer-typed `Column` of bin indices (`i8`..`u32`).
/// @param nbins  Number of bins; a non-negative integer `<= counts$length()`.
/// @param counts A numeric `Column` (length `>= nbins`) that receives the counts.
///   It is modified in place; the function returns `NULL`.
/// @return `NULL` (invisibly); the result is written into `counts`.
/// @examples
/// values <- allocate_scalar("u16", 6L)
/// values$set(c(0, 1, 1, 2, 2, 2))
/// counts <- allocate_scalar("i32", 3L)
/// bincount(values, 3L, counts)
/// counts$values()   # 1 2 3
/// @export
#[extendr]
fn bincount(values: &Column, nbins: i32, counts: &mut Column) {
    assert!(nbins >= 0, "`nbins` must be non-negative, got {nbins}");
    let nbins = nbins as usize;
    assert!(
        counts.len() >= nbins,
        "`counts` length ({}) must be at least `nbins` ({nbins})",
        counts.len()
    );

    // Recover the concrete element type and run the single generic histogram.
    let totals: Vec<u64> = match values.storage() {
        Storage::I8(v)  => histogram(v, nbins),
        Storage::U8(v)  => histogram(v, nbins),
        Storage::I16(v) => histogram(v, nbins),
        Storage::U16(v) => histogram(v, nbins),
        Storage::I32(v) => histogram(v, nbins),
        Storage::U32(v) => histogram(v, nbins),
        Storage::F32(_) | Storage::F64(_) =>
            panic!("`values` must be an integer-typed Column (bin indices), not float"),
    };

    // Zero the used bins of `counts`, then add the totals in.
    match counts.storage_mut() {
        Storage::I8(c)  => write_counts_u64(c, &totals, nbins),
        Storage::U8(c)  => write_counts_u64(c, &totals, nbins),
        Storage::I16(c) => write_counts_u64(c, &totals, nbins),
        Storage::U16(c) => write_counts_u64(c, &totals, nbins),
        Storage::I32(c) => write_counts_u64(c, &totals, nbins),
        Storage::U32(c) => write_counts_u64(c, &totals, nbins),
        Storage::F32(c) => write_counts_u64(c, &totals, nbins),
        Storage::F64(c) => write_counts_u64(c, &totals, nbins),
    }
}

// ── weights: any numeric value to add into a bin ────────────────────────────────

// `Sync` lets a shared `&[W]` be read from many Rayon threads. Weights may be
// signed, unsigned, or floating point — all widen to f64 for accumulation (exact
// for integer weights up to 2^53).
trait ToWeight: Copy + Sync {
    fn to_weight(self) -> f64;
}

macro_rules! impl_to_weight {
    ($($t:ty),*) => { $(
        impl ToWeight for $t {
            #[inline]
            fn to_weight(self) -> f64 { self as f64 }
        }
    )* };
}
impl_to_weight!(i8, u8, i16, u16, i32, u32, f32, f64);

// Weighted histogram: per Rayon task, accumulate `weights[i]` into bin
// `values[i]` in a private f64 buffer (no shared-bin collisions), then reduce the
// task buffers to a single per-bin total. `zip` pairs the two slices elementwise.
fn histogram_weighted<V: BinIndex, W: ToWeight>(
    values: &[V],
    weights: &[W],
    nbins: usize,
) -> Vec<f64> {
    values
        .par_iter()
        .zip(weights.par_iter())
        .fold(
            || vec![0.0f64; nbins],
            |mut local, (&v, &w)| {
                local[v.to_index()] += w.to_weight();
                local
            },
        )
        .reduce(
            || vec![0.0f64; nbins],
            |mut a, b| {
                for i in 0..nbins {
                    a[i] += b[i];
                }
                a
            },
        )
}

/// Weighted bincount: sum each element's weight into its bin (NumPy
/// `bincount(values, weights=...)`).
///
/// For each bin `b` in `0..nbins`, computes `sum of weights[i] over all i with
/// values[i] == b` and writes it into `counts[b]`. `values` must be an
/// integer-typed [Column] of bin indices in `0..nbins`; `weights` is any numeric
/// [Column] (signed, unsigned, or floating point) the SAME length as `values`;
/// `counts` is a numeric [Column] with at least `nbins` elements, **overwritten**
/// in its first `nbins` entries (entries at or beyond `nbins` are left untouched).
/// Parallelized with private per-thread accumulators, so there are no write
/// collisions.
///
/// @param values  An integer-typed `Column` of bin indices (`i8`..`u32`).
/// @param weights A numeric `Column` (any type), the same length as `values`.
/// @param nbins   Number of bins; a non-negative integer `<= counts$length()`.
/// @param counts  A numeric `Column` (length `>= nbins`) that receives the sums.
///   It is modified in place; the function returns `NULL`.
/// @return `NULL` (invisibly); the result is written into `counts`.
/// @examples
/// values  <- allocate_scalar("u16", 5L); values$set(c(0, 0, 1, 2, 2))
/// weights <- allocate_scalar("f64", 5L); weights$set(c(1.5, 2.5, 4, 1, 3))
/// counts  <- allocate_scalar("f64", 3L)
/// bincountw(values, weights, 3L, counts)
/// counts$values()   # 4 4 4
/// @export
#[extendr]
fn bincountw(values: &Column, weights: &Column, nbins: i32, counts: &mut Column) {
    assert!(nbins >= 0, "`nbins` must be non-negative, got {nbins}");
    let nbins = nbins as usize;
    assert!(
        counts.len() >= nbins,
        "`counts` length ({}) must be at least `nbins` ({nbins})",
        counts.len()
    );
    assert_eq!(
        values.len(), weights.len(),
        "`values` ({}) and `weights` ({}) must have the same length",
        values.len(), weights.len()
    );

    // Recover BOTH concrete element types (values' integer index type and weights'
    // numeric type) and call the single generic kernel. The two nested macros
    // expand, at compile time, into the value×weight match arms — so weights are
    // read in their native type with no copy, and only one `histogram_weighted` is
    // written. `dispatch_w` (inner) must be defined before `dispatch_vw` uses it.
    macro_rules! dispatch_w {
        ($vals:expr) => {
            match weights.storage() {
                Storage::I8(w)  => histogram_weighted($vals, w, nbins),
                Storage::U8(w)  => histogram_weighted($vals, w, nbins),
                Storage::I16(w) => histogram_weighted($vals, w, nbins),
                Storage::U16(w) => histogram_weighted($vals, w, nbins),
                Storage::I32(w) => histogram_weighted($vals, w, nbins),
                Storage::U32(w) => histogram_weighted($vals, w, nbins),
                Storage::F32(w) => histogram_weighted($vals, w, nbins),
                Storage::F64(w) => histogram_weighted($vals, w, nbins),
            }
        };
    }
    let totals: Vec<f64> = match values.storage() {
        Storage::I8(v)  => dispatch_w!(v),
        Storage::U8(v)  => dispatch_w!(v),
        Storage::I16(v) => dispatch_w!(v),
        Storage::U16(v) => dispatch_w!(v),
        Storage::I32(v) => dispatch_w!(v),
        Storage::U32(v) => dispatch_w!(v),
        Storage::F32(_) | Storage::F64(_) =>
            panic!("`values` must be an integer-typed Column (bin indices), not float"),
    };

    // Zero the used bins of `counts`, then add the weighted totals in.
    match counts.storage_mut() {
        Storage::I8(c)  => write_counts_f64(c, &totals, nbins),
        Storage::U8(c)  => write_counts_f64(c, &totals, nbins),
        Storage::I16(c) => write_counts_f64(c, &totals, nbins),
        Storage::U16(c) => write_counts_f64(c, &totals, nbins),
        Storage::I32(c) => write_counts_f64(c, &totals, nbins),
        Storage::U32(c) => write_counts_f64(c, &totals, nbins),
        Storage::F32(c) => write_counts_f64(c, &totals, nbins),
        Storage::F64(c) => write_counts_f64(c, &totals, nbins),
    }
}

extendr_module! {
    mod bincount;
    fn bincount;
    fn bincountw;
}
