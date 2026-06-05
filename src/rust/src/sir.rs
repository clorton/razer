// ════════════════════════════════════════════════════════════════════════════
// sir.rs — per-tick SIR step kernels operating on Column buffers.
//
// `sir_step` is the recovery half of an SIR tick: it advances infectious agents'
// countdown timers and moves the expired ones to Recovered, tallying recoveries
// per node into the current tick's slice of a report buffer. Later steps (e.g.
// transmission) will join it in the per-tick `run` loop.
//
// Orientation for readers coming from C / C++ / C# / Python:
//   * `state`, `timer`, `nodeid` are `Column` handles (Rust-owned, dtype-tagged
//     arrays); `as_u8_mut()` / `as_u16()` hand back a `&mut [u8]` / `&[u16]` view
//     straight into the buffer — mutating in place, no copy.
//   * `STATE_I as u8` casts the shared i32 state constant to the u8 the `state`
//     Column stores (state codes are small: S=0, I=2, R=3).
//   * `timer.saturating_sub(1)` subtracts but clamps at 0 instead of wrapping
//     (an unsigned underflow would otherwise jump to 255).
//   * `trait Increment` + the `impl_increment!` macro give one generic recovery
//     sweep that works for any numeric `recoveries` element type.
//
// This kernel is SERIAL on purpose: the per-agent state/timer writes are
// independent, but `recoveries[node] += 1` is a shared write that would race
// under naive parallelism. A parallel version would use per-thread node
// accumulators (cf. bincount.rs) — deferred until it's needed.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use crate::column::{Column, Storage};
use crate::epidemic::{STATE_I, STATE_R};

// Increment-by-one for a per-node recovery count of any numeric type.
trait Increment {
    fn inc(&mut self);
}

macro_rules! impl_increment {
    ($($t:ty),*) => { $(
        impl Increment for $t {
            #[inline]
            fn inc(&mut self) { *self += 1 as $t; }
        }
    )* };
}
impl_increment!(i8, u8, i16, u16, i32, u32, f32, f64);

// Recovery sweep over the first `count` agents. `rec` is the current tick's
// time-slice (length = node count); an agent that recovers bumps its node's entry.
fn recover<R: Increment>(
    state: &mut [u8],
    timer: &mut [u8],
    nodeid: &[u16],
    count: usize,
    rec: &mut [R],
) {
    let infectious = STATE_I as u8;
    let recovered = STATE_R as u8;
    for i in 0..count {
        if state[i] == infectious {
            timer[i] = timer[i].saturating_sub(1);
            if timer[i] == 0 {
                state[i] = recovered;
                rec[nodeid[i] as usize].inc();
            }
        }
    }
}

/// One SIR recovery step for a single tick.
///
/// For each of the first `count` agents whose `state` is infectious (`I`),
/// decrement its `timer` by one; when the timer reaches zero the agent becomes
/// recovered (`R`) and a recovery is tallied for its node in `recoveries` at this
/// `tick`. Susceptible and already-recovered agents are skipped.
///
/// `state` and `timer` are per-agent `u8` Columns and are mutated in place;
/// `nodeid` is the per-agent `u16` Column of 0-based node ids; `recoveries` is the
/// 2-D node report from [allocate_vector()] (shape `n_nodes × n_ticks`), of which
/// only the contiguous `n_nodes`-wide column for `tick` is written.
///
/// @param state      Per-agent `u8` state Column (mutated).
/// @param timer      Per-agent `u8` countdown Column (mutated).
/// @param nodeid     Per-agent `u16` 0-based node-id Column.
/// @param count      Number of active agents to process (e.g. `people$count`).
/// @param recoveries A 2-D node report Column (`n_nodes × n_ticks`); its `tick`
///   column receives the per-node recovery counts (mutated).
/// @param tick       0-based tick index selecting which column of `recoveries` to fill.
/// @return `NULL` (invisibly); `state`, `timer`, and `recoveries` are modified in place.
/// @export
#[extendr]
fn sir_step(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    recoveries: &mut Column,
    tick: i32,
) {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count = count as usize;
    let tick = tick as usize;

    // recoveries is a (n_ticks, n_nodes) report: each tick is a contiguous slice of
    // `slice_len` (= node count) elements; there are `n_slices` (= tick count) of them.
    let n_nodes = recoveries.slice_len();
    assert!(
        tick < recoveries.n_slices(),
        "`tick` ({tick}) out of range for recoveries with {} ticks",
        recoveries.n_slices()
    );
    // The current tick's contiguous slice within the slice-major buffer.
    let (start, end) = (tick * n_nodes, (tick + 1) * n_nodes);

    let state_s = state.as_u8_mut();
    let timer_s = timer.as_u8_mut();
    let node_s = nodeid.as_u16();
    assert!(
        count <= state_s.len() && count <= timer_s.len() && count <= node_s.len(),
        "`count` ({count}) exceeds the length of the people arrays"
    );

    // Dispatch on the recoveries element type; the generic sweep does the work.
    match recoveries.storage_mut() {
        Storage::I8(v)  => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
        Storage::U8(v)  => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
        Storage::I16(v) => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
        Storage::U16(v) => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
        Storage::I32(v) => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
        Storage::U32(v) => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
        Storage::F32(v) => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
        Storage::F64(v) => recover(state_s, timer_s, node_s, count, &mut v[start..end]),
    }
}

extendr_module! {
    mod sir;
    fn sir_step;
}
