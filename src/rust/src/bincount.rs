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
/// `counts` is a numeric [Column] and the result is written into **slice `slot`**
/// of it — the `nbins` entries `slot*slice_len .. slot*slice_len + nbins` —
/// overwriting them while leaving the rest of the slice (and other slices)
/// untouched. For a scalar `counts` (shape `(1, n)`) the only slice is `slot = 0`,
/// the whole vector; for a 2-D report (e.g. `(n_ticks, n_nodes)`) `slot` selects a
/// tick's row. The work is parallelized with private per-thread histograms, so
/// there are no write collisions.
///
/// @param values An integer-typed `Column` of bin indices (`i8`..`u32`).
/// @param nbins  Number of bins; a non-negative integer `<= counts`'s slice length.
/// @param counts A numeric `Column` that receives the counts (modified in place).
/// @param slot   Which slice of `counts` to write; a non-negative integer
///   `< counts`'s slice count. Defaults to `0` (the whole vector for a scalar
///   `counts`, or the first tick of a report). `@param` is optional.
/// @return `NULL` (invisibly); the result is written into `counts`.
/// @examples
/// values <- allocate_scalar("u16", 6L)
/// values$set(c(0, 1, 1, 2, 2, 2))
/// counts <- allocate_scalar("i32", 3L)
/// bincount(values, 3L, counts)
/// counts$values()   # 1 2 3
/// @noRd
#[extendr]
fn bincount_impl(values: &Column, nbins: i32, counts: &mut Column, slot: i32) {
    assert!(nbins >= 0, "`nbins` must be non-negative, got {nbins}");
    let nbins = nbins as usize;
    assert!(slot >= 0, "`slot` must be non-negative, got {slot}");
    let slot = slot as usize;
    assert!(
        slot < counts.n_slices(),
        "`slot` ({slot}) out of range for counts with {} slices",
        counts.n_slices()
    );
    let slice_len = counts.slice_len();
    assert!(
        nbins <= slice_len,
        "`counts` slice length ({slice_len}) must be at least `nbins` ({nbins})"
    );
    let start = slot * slice_len;

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

    // Zero slice `slot`'s first `nbins` entries, then add the totals in. The
    // sub-slice `[start .. start+nbins]` is the destination within `counts`.
    let end = start + nbins;
    match counts.storage_mut() {
        Storage::I8(c)  => write_counts_u64(&mut c[start..end], &totals, nbins),
        Storage::U8(c)  => write_counts_u64(&mut c[start..end], &totals, nbins),
        Storage::I16(c) => write_counts_u64(&mut c[start..end], &totals, nbins),
        Storage::U16(c) => write_counts_u64(&mut c[start..end], &totals, nbins),
        Storage::I32(c) => write_counts_u64(&mut c[start..end], &totals, nbins),
        Storage::U32(c) => write_counts_u64(&mut c[start..end], &totals, nbins),
        Storage::F32(c) => write_counts_u64(&mut c[start..end], &totals, nbins),
        Storage::F64(c) => write_counts_u64(&mut c[start..end], &totals, nbins),
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
/// values[i] == b` and writes it into **slice `slot`** of `counts` (the entries
/// `slot*slice_len .. slot*slice_len + nbins`), overwriting them. `values` must be
/// an integer-typed [Column] of bin indices in `0..nbins`; `weights` is any numeric
/// [Column] (signed, unsigned, or floating point) the SAME length as `values`;
/// `counts` is a numeric [Column] whose slice length is `>= nbins`. For a scalar
/// `counts` the only slice is `slot = 0` (the whole vector); for a 2-D report
/// `slot` selects a tick's row. Parallelized with private per-thread accumulators,
/// so there are no write collisions.
///
/// @param values  An integer-typed `Column` of bin indices (`i8`..`u32`).
/// @param weights A numeric `Column` (any type), the same length as `values`.
/// @param nbins   Number of bins; a non-negative integer `<= counts`'s slice length.
/// @param counts  A numeric `Column` that receives the sums (modified in place).
/// @param slot    Which slice of `counts` to write; a non-negative integer
///   `< counts`'s slice count. Defaults to `0`.
/// @return `NULL` (invisibly); the result is written into `counts`.
/// @examples
/// values  <- allocate_scalar("u16", 5L); values$set(c(0, 0, 1, 2, 2))
/// weights <- allocate_scalar("f64", 5L); weights$set(c(1.5, 2.5, 4, 1, 3))
/// counts  <- allocate_scalar("f64", 3L)
/// bincountw(values, weights, 3L, counts)
/// counts$values()   # 4 4 4
/// @noRd
#[extendr]
fn bincountw_impl(values: &Column, weights: &Column, nbins: i32, counts: &mut Column, slot: i32) {
    assert!(nbins >= 0, "`nbins` must be non-negative, got {nbins}");
    let nbins = nbins as usize;
    assert!(slot >= 0, "`slot` must be non-negative, got {slot}");
    let slot = slot as usize;
    assert!(
        slot < counts.n_slices(),
        "`slot` ({slot}) out of range for counts with {} slices",
        counts.n_slices()
    );
    let slice_len = counts.slice_len();
    assert!(
        nbins <= slice_len,
        "`counts` slice length ({slice_len}) must be at least `nbins` ({nbins})"
    );
    assert_eq!(
        values.len(), weights.len(),
        "`values` ({}) and `weights` ({}) must have the same length",
        values.len(), weights.len()
    );
    let start = slot * slice_len;

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

    // Zero slice `slot`'s first `nbins` entries, then add the weighted totals in.
    let end = start + nbins;
    match counts.storage_mut() {
        Storage::I8(c)  => write_counts_f64(&mut c[start..end], &totals, nbins),
        Storage::U8(c)  => write_counts_f64(&mut c[start..end], &totals, nbins),
        Storage::I16(c) => write_counts_f64(&mut c[start..end], &totals, nbins),
        Storage::U16(c) => write_counts_f64(&mut c[start..end], &totals, nbins),
        Storage::I32(c) => write_counts_f64(&mut c[start..end], &totals, nbins),
        Storage::U32(c) => write_counts_f64(&mut c[start..end], &totals, nbins),
        Storage::F32(c) => write_counts_f64(&mut c[start..end], &totals, nbins),
        Storage::F64(c) => write_counts_f64(&mut c[start..end], &totals, nbins),
    }
}

// ── count_by_where: predicate-filtered, count-aware bincount ─────────────────────
//
// "How many agents in each node are in state E?" or "...are under five (tick - dob <
// 5*365)?" Both are a bincount of `group` (e.g. nodeid) restricted to the agents whose
// some-property `prop` satisfies a comparison against a threshold `value`, scanning only
// the first `count` agents (the live prefix — capacity-sized Columns may trail inactive
// slots). One parallel pass, no intermediate mask Column, no copy of `prop` to R.

// Build a per-group histogram counting only elements where `pred(prop[i])` holds.
// `pred` is a non-capturing comparison baked from the op + threshold, so it is `Copy`
// and trivially `Send + Sync`. Same private-buffer-then-reduce trick as `histogram`.
fn histogram_filtered<V: BinIndex, P: ToWeight, F: Fn(f64) -> bool + Sync + Send>(
    group: &[V],
    prop: &[P],
    nbins: usize,
    pred: F,
) -> Vec<u64> {
    group
        .par_iter()
        .zip(prop.par_iter())
        .fold(
            || vec![0u64; nbins],
            |mut local, (&g, &p)| {
                if pred(p.to_weight()) {
                    local[g.to_index()] += 1;
                }
                local
            },
        )
        .reduce(
            || vec![0u64; nbins],
            |mut a, b| {
                for i in 0..nbins {
                    a[i] += b[i];
                }
                a
            },
        )
}

/// Count, per group, the agents whose property satisfies a comparison (filtered bincount).
///
/// For each group `g` in `0..n_groups`, counts how many of the first `count` agents both
/// have `group[i] == g` AND satisfy `prop[i] <op> value`, and writes the totals into
/// **slice `slot`** of `counts` (entries `slot*slice_len .. slot*slice_len + n_groups`),
/// overwriting them. This is `bincount` restricted to a predicate and to the live prefix
/// `0..count` — e.g. with `group = nodeid`, `prop = state`, `op = "eq"`, `value = E` it
/// tallies the exposed by node; with `prop = dob`, `op = "gt"`, `value = tick - 5*365` it
/// tallies the under-fives by node (since `dob = -age`). `group` must be an integer-typed
/// [Column] of indices in `0..n_groups`; `prop` is any numeric [Column] (compared as
/// `f64`); `counts` is numeric with slice length `>= n_groups`. Parallelized with private
/// per-thread histograms, so there are no write collisions.
///
/// @param group    An integer-typed `Column` of group indices (`i8`..`u32`), e.g. nodeid.
/// @param n_groups Number of groups; a non-negative integer `<= counts`'s slice length.
/// @param prop     A numeric `Column` (any type) holding the per-agent property to test.
/// @param op       Comparison: one of `"eq"`, `"ne"`, `"lt"`, `"le"`, `"gt"`, `"ge"`.
/// @param value    The threshold the property is compared against (a double).
/// @param count    How many leading agents to scan (the active count); `<= group`/`prop`
///   length.
/// @param counts   A numeric `Column` that receives the per-group totals (modified in place).
/// @param slot     Which slice of `counts` to write; a non-negative integer
///   `< counts`'s slice count. Defaults to `0`.
/// @return `NULL` (invisibly); the result is written into `counts`.
/// @noRd
#[extendr]
fn count_by_where_impl(
    group: &Column,
    n_groups: i32,
    prop: &Column,
    op: &str,
    value: f64,
    count: i32,
    counts: &mut Column,
    slot: i32,
) {
    assert!(n_groups >= 0, "`n_groups` must be non-negative, got {n_groups}");
    let n_groups = n_groups as usize;
    assert!(slot >= 0, "`slot` must be non-negative, got {slot}");
    let slot = slot as usize;
    assert!(
        slot < counts.n_slices(),
        "`slot` ({slot}) out of range for counts with {} slices",
        counts.n_slices()
    );
    let slice_len = counts.slice_len();
    assert!(
        n_groups <= slice_len,
        "`counts` slice length ({slice_len}) must be at least `n_groups` ({n_groups})"
    );
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    let count = count as usize;
    assert!(
        count <= group.len() && count <= prop.len(),
        "`count` ({count}) exceeds `group` ({}) or `prop` ({}) length",
        group.len(), prop.len()
    );

    // Bake the comparison into a non-capturing fn pointer (all arms share one type), then
    // wrap it with the captured threshold. `move` copies `value` into the closure.
    let cmp: fn(f64, f64) -> bool = match op {
        "eq" => |a, b| a == b,
        "ne" => |a, b| a != b,
        "lt" => |a, b| a < b,
        "le" => |a, b| a <= b,
        "gt" => |a, b| a > b,
        "ge" => |a, b| a >= b,
        _ => panic!("`op` must be one of \"eq\", \"ne\", \"lt\", \"le\", \"gt\", \"ge\", got {op:?}"),
    };
    let pred = move |x: f64| cmp(x, value);

    // Recover both element types (group's integer index type, prop's numeric type) over
    // the live prefix `[..count]` and run the single generic filtered histogram. `pred`
    // is `Copy`, so each arm gets its own copy.
    macro_rules! dispatch_prop {
        ($g:expr) => {
            match prop.storage() {
                Storage::I8(p)  => histogram_filtered($g, &p[..count], n_groups, pred),
                Storage::U8(p)  => histogram_filtered($g, &p[..count], n_groups, pred),
                Storage::I16(p) => histogram_filtered($g, &p[..count], n_groups, pred),
                Storage::U16(p) => histogram_filtered($g, &p[..count], n_groups, pred),
                Storage::I32(p) => histogram_filtered($g, &p[..count], n_groups, pred),
                Storage::U32(p) => histogram_filtered($g, &p[..count], n_groups, pred),
                Storage::F32(p) => histogram_filtered($g, &p[..count], n_groups, pred),
                Storage::F64(p) => histogram_filtered($g, &p[..count], n_groups, pred),
            }
        };
    }
    let totals: Vec<u64> = match group.storage() {
        Storage::I8(g)  => dispatch_prop!(&g[..count]),
        Storage::U8(g)  => dispatch_prop!(&g[..count]),
        Storage::I16(g) => dispatch_prop!(&g[..count]),
        Storage::U16(g) => dispatch_prop!(&g[..count]),
        Storage::I32(g) => dispatch_prop!(&g[..count]),
        Storage::U32(g) => dispatch_prop!(&g[..count]),
        Storage::F32(_) | Storage::F64(_) =>
            panic!("`group` must be an integer-typed Column (group indices), not float"),
    };

    let start = slot * slice_len;
    let end = start + n_groups;
    match counts.storage_mut() {
        Storage::I8(c)  => write_counts_u64(&mut c[start..end], &totals, n_groups),
        Storage::U8(c)  => write_counts_u64(&mut c[start..end], &totals, n_groups),
        Storage::I16(c) => write_counts_u64(&mut c[start..end], &totals, n_groups),
        Storage::U16(c) => write_counts_u64(&mut c[start..end], &totals, n_groups),
        Storage::I32(c) => write_counts_u64(&mut c[start..end], &totals, n_groups),
        Storage::U32(c) => write_counts_u64(&mut c[start..end], &totals, n_groups),
        Storage::F32(c) => write_counts_u64(&mut c[start..end], &totals, n_groups),
        Storage::F64(c) => write_counts_u64(&mut c[start..end], &totals, n_groups),
    }
}

extendr_module! {
    mod bincount;
    fn bincount_impl;
    fn bincountw_impl;
    fn count_by_where_impl;
}
