// ════════════════════════════════════════════════════════════════════════════
// measles.rs — the combined timed-transition kernel for the measles model.
//
// Measles tracks three TIMED compartments, each counted down by the same per-agent
// `timer`: maternal immunity (M), the exposed/latent period (E), and the infectious
// period (I). `measles_step` advances all three in a SINGLE pass over the agents:
//
//   M -> S   maternal antibodies wane (timer expires) -> Susceptible (untimed).
//   E -> I   incubation ends -> Infectious; a fresh infectious-period timer is drawn.
//   I -> R   infectiousness ends -> Recovered (untimed; measles confers lifelong immunity).
//
// Doing all three in one pass — branching on each agent's state at the START of the
// pass — means every agent is touched exactly once per tick, so an agent that does
// E->I this tick is NOT then processed as I->R in the same tick. That is the
// downstream-first ordering CLAUDE.md describes, achieved structurally rather than by
// sequencing separate kernels.
//
// The timer is a uint16 because maternal immunity lasts ~270 days, beyond a uint8's 255
// range. Decrement is GUARDED (`if t > 0 { t -= 1 }`) so a u16 never underflows.
//
// Parallel across agents with a private per-task per-node tally reduced at the end (see
// the project memory). RNG (for the E->I infectious-period draw) is thread-local.
// State codes live in epidemic.rs; the u8 `state` Column stores D (= -1) as 255.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::column::Column;
use crate::distributions::Distribution;
use crate::epidemic::{STATE_S, STATE_E, STATE_I, STATE_R, STATE_M};

// Round a distribution draw to a whole-tick u16 timer, clamped to [1, u16::MAX]: a
// freshly entered timed state lasts at least one tick.
#[inline]
fn timer_u16(x: f64) -> u16 {
    x.round().max(1.0).min(u16::MAX as f64) as u16
}

/// Advance the measles timed compartments for one tick (M→S, E→I, I→R).
///
/// For each of the first `count` agents, decrements its `timer` (uint16) and, on expiry,
/// transitions it: M→S (timer left at 0), E→I (timer set to a fresh draw from
/// `inf_duration`, the infectious period), or I→R (timer left at 0). Each agent is in at
/// most one timed state, so it is processed at most once. The per-compartment census
/// deltas are applied at column `tick + 1` (the working columns the caller has already
/// carried forward): M and E lose their leavers, R gains its arrivals, S gains the M→S
/// waners, and I gains the E→I arrivals minus the I→R leavers. Parallelized across cores
/// with private per-thread node buffers reduced at the end.
///
/// @param state   Per-agent `u8` state Column (mutated).
/// @param timer   Per-agent `u16` countdown Column (mutated; decremented, reset on E→I).
/// @param nodeid  Per-agent `u16` 0-based node-id Column.
/// @param count   Number of active agents to process.
/// @param m_count,s_count,e_count,i_count,r_count  `n_ticks x n_nodes` i32 census
///   Columns kept in sync (mutated at column `tick + 1`).
/// @param inf_duration A `Distribution` for the infectious period drawn on E→I.
/// @param tick    0-based tick index.
/// @return `NULL` (invisibly); the Columns are modified in place.
/// @export
#[extendr]
fn measles_step(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    m_count: &mut Column,
    s_count: &mut Column,
    e_count: &mut Column,
    i_count: &mut Column,
    r_count: &mut Column,
    inf_duration: &Distribution,
    tick: i32,
) {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count = count as usize;
    let tick  = tick as usize;

    let n = m_count.slice_len();
    for (nm, c) in [("m", &*m_count), ("s", &*s_count), ("e", &*e_count),
                    ("i", &*i_count), ("r", &*r_count)] {
        assert_eq!(c.slice_len(), n, "`{nm}_count` slice length must equal n_nodes");
    }
    assert!(tick + 1 < m_count.n_slices(), "`tick`+1 out of range for the census buffers");

    let m_code = STATE_M as u8;
    let s_code = STATE_S as u8;
    let e_code = STATE_E as u8;
    let i_code = STATE_I as u8;
    let r_code = STATE_R as u8;

    let nthreads = rayon::current_num_threads().max(1);
    let chunk = ((count + nthreads - 1) / nthreads).max(1);
    let st  = &mut state.as_u8_mut()[..count];
    let tm  = &mut timer.as_u16_mut()[..count];
    let nid = &nodeid.as_u16()[..count];

    // Per-node tally of transitions: slot 0 = M->S, 1 = E->I, 2 = I->R.
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|((s_chunk, t_chunk), node_chunk)| {
            let mut rng = rand::thread_rng();    // per-worker RNG for the E->I draw
            let mut local = vec![0i64; 3 * n];
            for j in 0..s_chunk.len() {
                let k = node_chunk[j] as usize;
                if s_chunk[j] == m_code {
                    if t_chunk[j] > 0 { t_chunk[j] -= 1; }      // guarded: no u16 underflow
                    if t_chunk[j] == 0 {
                        s_chunk[j] = s_code;                    // M -> S (untimed)
                        local[k] += 1;
                    }
                } else if s_chunk[j] == e_code {
                    if t_chunk[j] > 0 { t_chunk[j] -= 1; }
                    if t_chunk[j] == 0 {
                        s_chunk[j] = i_code;                    // E -> I
                        t_chunk[j] = timer_u16(inf_duration.sample(&mut rng));
                        local[n + k] += 1;
                    }
                } else if s_chunk[j] == i_code {
                    if t_chunk[j] > 0 { t_chunk[j] -= 1; }
                    if t_chunk[j] == 0 {
                        s_chunk[j] = r_code;                    // I -> R (untimed)
                        local[2 * n + k] += 1;
                    }
                }
                // S, R, D: no timer, nothing to do.
            }
            local
        })
        .reduce(|| vec![0i64; 3 * n], |mut a, b| {
            for k in 0..3 * n { a[k] += b[k]; }
            a
        });

    // Apply the census deltas at column tick+1.
    let dst = (tick + 1) * n;
    { let m = m_count.as_i32_mut(); for k in 0..n { m[dst + k] -= tally[k]         as i32; } } // M loses waners
    { let s = s_count.as_i32_mut(); for k in 0..n { s[dst + k] += tally[k]         as i32; } } // S gains waners
    { let e = e_count.as_i32_mut(); for k in 0..n { e[dst + k] -= tally[n + k]     as i32; } } // E loses onsets
    // I gains the E->I onsets and loses the I->R recoveries.
    { let i = i_count.as_i32_mut(); for k in 0..n { i[dst + k] += (tally[n + k] - tally[2 * n + k]) as i32; } }
    { let r = r_count.as_i32_mut(); for k in 0..n { r[dst + k] += tally[2 * n + k] as i32; } } // R gains recoveries
}

/// Transmission step (S→`to_state`) writing a uint16 timer — the measles counterpart of
/// [transmission()], which writes a u8 timer.
///
/// Identical to [transmission()] in every respect except the `timer` Column is `u16`
/// (clamped to `[1, 65535]`): converts column `tick` of `foi` to a per-node probability
/// `1 - exp(-foi)` once per node, then for each susceptible agent moves it (with that
/// probability) into `to_state` — `E` for measles — drawing its `timer` from `duration`
/// (the incubation period). New cases per node are applied as a delta to column `tick+1`
/// of the `s_count`/`to_count` census and recorded in column `tick` of `incidence`.
///
/// @param state      Per-agent `u8` state Column (mutated).
/// @param timer      Per-agent `u16` countdown Column for the receiving state (mutated).
/// @param nodeid     Per-agent `u16` 0-based node-id Column.
/// @param count      Number of active agents to process.
/// @param foi        `n_ticks x n_nodes` f64 FOI Column (from [calc_foi()]); column `tick` read.
/// @param s_count    i32 census Column for `S` (mutated at `tick+1`).
/// @param to_count   i32 census Column for the receiving state (mutated at `tick+1`).
/// @param incidence  i32 flow Column (mutated at `tick`).
/// @param tick       0-based tick index.
/// @param to_state   State code new cases enter (e.g. `laser_states()[["E"]]`).
/// @param duration   A Distribution from which the receiving state's timer is drawn.
/// @return `NULL` (invisibly); the Columns are modified in place.
/// @export
#[extendr]
fn transmission_u16(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    foi: &Column,
    s_count: &mut Column,
    to_count: &mut Column,
    incidence: &mut Column,
    tick: i32,
    to_state: i32,
    duration: &Distribution,
) {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    assert!((0..=255).contains(&to_state), "`to_state` must be a state code in [0, 255], got {to_state}");
    let count = count as usize;
    let tick = tick as usize;
    let to = to_state as u8;

    let n = foi.slice_len();
    assert_eq!(s_count.slice_len(), n, "census Columns must share n_nodes");
    assert_eq!(to_count.slice_len(), n, "census Columns must share n_nodes");
    assert_eq!(incidence.slice_len(), n, "`incidence` slice length must equal n_nodes");
    assert!(tick < foi.n_slices(), "`tick` out of range for `foi`");
    assert!(tick + 1 < to_count.n_slices(), "`tick`+1 out of range for the census buffers");
    assert!(tick < incidence.n_slices(), "`tick` out of range for `incidence`");

    // Per-NODE infection probability (one exp() per node, reused for every agent there).
    let foi_col = &foi.as_f64()[tick * n..(tick + 1) * n];
    let p: Vec<f64> = foi_col.iter().map(|&lambda| 1.0 - (-lambda).exp()).collect();

    let susceptible = STATE_S as u8;
    let nthreads = rayon::current_num_threads().max(1);
    let chunk = ((count + nthreads - 1) / nthreads).max(1);
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u16_mut()[..count];   // u16 timer (vs. transmission()'s u8)
    let nid = &nodeid.as_u16()[..count];
    let new: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|((s, t), node)| {
            let mut rng = rand::thread_rng();
            let mut local = vec![0i64; n];
            for j in 0..s.len() {
                if s[j] == susceptible {
                    let k = node[j] as usize;
                    if p[k] > 0.0 && rng.gen::<f64>() < p[k] {
                        s[j] = to;
                        t[j] = timer_u16(duration.sample(&mut rng));
                        local[k] += 1;
                    }
                }
            }
            local
        })
        .reduce(|| vec![0i64; n], |mut a, b| {
            for k in 0..n { a[k] += b[k]; }
            a
        });

    let dst = (tick + 1) * n;
    { let sc = s_count.as_i32_mut();  for k in 0..n { sc[dst + k] -= new[k] as i32; } }
    { let tc = to_count.as_i32_mut(); for k in 0..n { tc[dst + k] += new[k] as i32; } }
    { let flow = incidence.as_i32_mut(); for k in 0..n { flow[tick * n + k] = new[k] as i32; } }
}

extendr_module! {
    mod measles;
    fn measles_step;
    fn transmission_u16;
}
