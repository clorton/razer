// ════════════════════════════════════════════════════════════════════════════
// rng.rs — reproducible, seedable RNG for the kernels.
//
// All kernel randomness flows through here so that a single `set_seed()` from R makes
// an entire razer run reproducible — independent of the machine's CPU/thread count.
//
// How reproducibility survives Rayon parallelism: the parallel kernels split the agents
// into FIXED-size chunks (`RNG_CHUNK`, NOT one-per-thread), and chunk `ci` is given an
// RNG seeded deterministically from (call base, ci). So chunk `ci` always covers the same
// agent range and draws the same stream no matter how Rayon schedules chunks onto
// threads. Per-node count reductions are integer sums (order-independent) and each agent
// is mutated by exactly one chunk, so the whole result is a pure function of the seed.
//
// Each kernel invocation pulls ONE `next_call_base()` on the (single-threaded) R calling
// thread before fanning out; when seeded, that advances a deterministic counter so
// successive kernel calls and ticks get independent streams. When NOT seeded (no
// `set_seed`), the base is drawn from OS entropy, so the default behaviour is a fresh
// random run — exactly as before.
//
// Orientation for readers coming from C / C++ / C# / Python: `AtomicU64`/`AtomicBool`
// are lock-free shared globals (like `std::atomic` / `Interlocked`); `SmallRng` is rand's
// fast non-cryptographic generator (xoshiro256++); `seed_from_u64` deterministically
// seeds it. `#[extendr]` exposes `set_seed`/`unset_seed` to R.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rand::{SeedableRng, RngCore};
use rand::rngs::SmallRng;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

/// The model RNG: fast, non-cryptographic, seedable.
pub(crate) type ModelRng = SmallRng;

/// Fixed per-chunk agent count for the parallel kernels (independent of thread count, so
/// chunking — and therefore the seeded streams — is reproducible across machines).
pub(crate) const RNG_CHUNK: usize = 16_384;

static SEED: AtomicU64 = AtomicU64::new(0);
static COUNTER: AtomicU64 = AtomicU64::new(0);
static SEEDED: AtomicBool = AtomicBool::new(false);

const GOLDEN: u64 = 0x9e37_79b9_7f4a_7c15;

// splitmix64 finalizer — decorrelates nearby seeds.
#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// A per-kernel-call base seed. Call ONCE per kernel invocation, on the single-threaded
/// R calling thread, before the parallel fan-out. Seeded: deterministic counter advance.
/// Unseeded: fresh entropy (random run).
pub(crate) fn next_call_base() -> u64 {
    if SEEDED.load(Ordering::Acquire) {
        let s = SEED.load(Ordering::Relaxed);
        let c = COUNTER.fetch_add(1, Ordering::Relaxed);
        mix64(s ^ c.wrapping_mul(GOLDEN))
    } else {
        SmallRng::from_entropy().next_u64()
    }
}

/// A per-chunk RNG seeded deterministically from the call base and chunk index.
#[inline]
pub(crate) fn chunk_rng(base: u64, chunk_index: usize) -> ModelRng {
    SmallRng::seed_from_u64(mix64(base ^ (chunk_index as u64).wrapping_mul(GOLDEN)))
}

/// A single RNG for sequential kernels (one call base, chunk 0).
pub(crate) fn single_rng() -> ModelRng {
    chunk_rng(next_call_base(), 0)
}

/// Set the global random seed, making subsequent razer runs reproducible.
///
/// After `set_seed(s)`, every kernel's randomness is a deterministic function of `s` and
/// the order of kernel calls — identical on every run and every machine, regardless of
/// CPU/thread count. Call it once at the start of a script. `unset_seed()` reverts to a
/// fresh (entropy-seeded) RNG.
///
/// @param seed A finite, non-negative number (used as a 64-bit seed).
/// @return `NULL`, invisibly.
/// @examples
/// set_seed(42)
/// @export
#[extendr]
fn set_seed(seed: f64) {
    assert!(seed.is_finite() && seed >= 0.0, "`seed` must be a finite, non-negative number, got {seed}");
    SEED.store(seed as u64, Ordering::Relaxed);
    COUNTER.store(0, Ordering::Relaxed);
    SEEDED.store(true, Ordering::Release);
}

/// Revert to a non-reproducible, entropy-seeded RNG (undo [set_seed()]).
///
/// @return `NULL`, invisibly.
/// @export
#[extendr]
fn unset_seed() {
    SEEDED.store(false, Ordering::Release);
}

extendr_module! {
    mod rng;
    fn set_seed;
    fn unset_seed;
}
