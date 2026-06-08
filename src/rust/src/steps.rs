// ════════════════════════════════════════════════════════════════════════════
// steps.rs — the timed-transition step kernels for the model menagerie.
//
// Three kernels cover all eight SI/SEI/SIS/SEIS/SIR/SEIR/SIRS/SEIRS models. Each is a
// SINGLE pass over the agents, branching on each agent's state at the start of the pass
// (so a just-arrived timed state is never decremented again the same tick — the
// downstream-first ordering, achieved structurally). Each leads with M→S maternal
// waning, so a maternal-immunity compartment can be added to any model.
//
//   step_si   : M→S, E→I                      → SI, SEI
//   step_sir  : M→S, E→I, I→{S|R} (param)      → SIS, SIR, SEIS, SEIR
//   step_sirs : M→S, E→I, I→R(+imm), R→S       → SIRS, SEIRS
//
// Following the project convention, the kernels touch NO node census/flow buffers: they
// mutate the per-agent `state`/`timer` arrays and RETURN per-node transition counts as a
// named list of integer vectors. The model applies the counts to whichever compartments
// it maintains (a model with no E ignores `onset`, never allocating an E census). Timers
// are u16 (maternal/immunity periods exceed a u8's 255). Parallel across agents (Rayon)
// with private per-node tallies reduced at the end; RNG is thread-local.
//
// Returned list elements (each an integer[n_nodes]): `waned` = M→S, `onset` = E→I,
// `cleared` = I→absorbing (step_sir), `recovered` = I→R, `waned_r` = R→S (step_sirs).
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use crate::column::Column;
use crate::distributions::Distribution;
use crate::epidemic::{STATE_S, STATE_E, STATE_I, STATE_R, STATE_M};
use crate::rng;

#[inline]
fn timer_u16(x: f64) -> u16 {
    x.round().max(1.0).min(u16::MAX as f64) as u16
}

// Convert one per-node block of the tally to an R integer vector.
fn block(tally: &[i64], slot: usize, n: usize) -> Vec<i32> {
    tally[slot * n..(slot + 1) * n].iter().map(|&x| x as i32).collect()
}

/// Advance M→S (maternal waning) and E→I (incubation) for one tick — SI / SEI.
///
/// `I` is terminal (no exit). Returns `list(waned, onset)` of per-node counts.
///
/// @param state   Per-agent `u8` state Column (mutated).
/// @param timer   Per-agent `u16` timer Column (mutated; E→I draws an infectious timer).
/// @param nodeid  Per-agent `u16` 0-based node-id Column.
/// @param count   Number of active agents to process.
/// @param n_nodes Number of nodes (the length of each returned vector).
/// @param inf_duration A Distribution for the infectious period set on E→I.
/// @return `list(waned = integer[n_nodes], onset = integer[n_nodes])`.
/// @export
#[extendr]
fn step_si(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    n_nodes: i32,
    inf_duration: &Distribution,
) -> List {
    assert!(count >= 0 && n_nodes >= 0, "`count`/`n_nodes` must be non-negative");
    let count = count as usize;
    let n = n_nodes as usize;
    let (m_code, s_code, e_code, i_code) =
        (STATE_M as u8, STATE_S as u8, STATE_E as u8, STATE_I as u8);

    let base = rng::next_call_base();
    let chunk = rng::RNG_CHUNK;
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u16_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .enumerate()
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|(((ci, s), t), node)| {
            let mut rng = rng::chunk_rng(base, ci);
            let mut local = vec![0i64; 2 * n];
            for j in 0..s.len() {
                let k = node[j] as usize;
                if s[j] == m_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 { s[j] = s_code; local[k] += 1; }
                } else if s[j] == e_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 {
                        s[j] = i_code;
                        t[j] = timer_u16(inf_duration.sample(&mut rng));
                        local[n + k] += 1;
                    }
                }
            }
            local
        })
        .reduce(|| vec![0i64; 2 * n], |mut a, b| { for k in 0..2 * n { a[k] += b[k]; } a });
    list!(waned = block(&tally, 0, n), onset = block(&tally, 1, n))
}

/// Advance M→S, E→I, and I→`absorbing_state` for one tick — SIS / SIR / SEIS / SEIR.
///
/// `absorbing_state` is the untimed destination of the infectious period: `S` (SIS/SEIS)
/// or `R` (SIR/SEIR). Returns `list(waned, onset, cleared)` of per-node counts (`cleared`
/// is the I→`absorbing_state` flow).
///
/// @param state,timer,nodeid,count,n_nodes  As in [step_si()].
/// @param inf_duration   A Distribution for the infectious period set on E→I.
/// @param absorbing_state State code I clears to (`laser_states()[["S"]]` or `[["R"]]`).
/// @return `list(waned, onset, cleared)` of `integer[n_nodes]`.
/// @export
#[extendr]
fn step_sir(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    n_nodes: i32,
    inf_duration: &Distribution,
    absorbing_state: i32,
) -> List {
    assert!(count >= 0 && n_nodes >= 0, "`count`/`n_nodes` must be non-negative");
    assert!((0..=255).contains(&absorbing_state), "`absorbing_state` must be in [0, 255]");
    let count = count as usize;
    let n = n_nodes as usize;
    let (m_code, s_code, e_code, i_code) =
        (STATE_M as u8, STATE_S as u8, STATE_E as u8, STATE_I as u8);
    let abs = absorbing_state as u8;

    let base = rng::next_call_base();
    let chunk = rng::RNG_CHUNK;
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u16_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .enumerate()
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|(((ci, s), t), node)| {
            let mut rng = rng::chunk_rng(base, ci);
            let mut local = vec![0i64; 3 * n];
            for j in 0..s.len() {
                let k = node[j] as usize;
                if s[j] == m_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 { s[j] = s_code; local[k] += 1; }
                } else if s[j] == e_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 {
                        s[j] = i_code;
                        t[j] = timer_u16(inf_duration.sample(&mut rng));
                        local[n + k] += 1;
                    }
                } else if s[j] == i_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 { s[j] = abs; local[2 * n + k] += 1; }   // untimed destination
                }
            }
            local
        })
        .reduce(|| vec![0i64; 3 * n], |mut a, b| { for k in 0..3 * n { a[k] += b[k]; } a });
    list!(waned = block(&tally, 0, n), onset = block(&tally, 1, n), cleared = block(&tally, 2, n))
}

/// Advance M→S, E→I, I→R (with waning immunity), and R→S for one tick — SIRS / SEIRS.
///
/// On I→R a fresh immunity timer is drawn from `imm_duration`; R→S fires when it expires.
/// Returns `list(waned, onset, recovered, waned_r)` of per-node counts.
///
/// @param state,timer,nodeid,count,n_nodes  As in [step_si()].
/// @param inf_duration A Distribution for the infectious period set on E→I.
/// @param imm_duration A Distribution for the immunity period set on I→R.
/// @return `list(waned, onset, recovered, waned_r)` of `integer[n_nodes]`.
/// @export
#[extendr]
fn step_sirs(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    n_nodes: i32,
    inf_duration: &Distribution,
    imm_duration: &Distribution,
) -> List {
    assert!(count >= 0 && n_nodes >= 0, "`count`/`n_nodes` must be non-negative");
    let count = count as usize;
    let n = n_nodes as usize;
    let (m_code, s_code, e_code, i_code, r_code) =
        (STATE_M as u8, STATE_S as u8, STATE_E as u8, STATE_I as u8, STATE_R as u8);

    let base = rng::next_call_base();
    let chunk = rng::RNG_CHUNK;
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u16_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .enumerate()
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|(((ci, s), t), node)| {
            let mut rng = rng::chunk_rng(base, ci);
            let mut local = vec![0i64; 4 * n];
            for j in 0..s.len() {
                let k = node[j] as usize;
                if s[j] == m_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 { s[j] = s_code; local[k] += 1; }
                } else if s[j] == e_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 {
                        s[j] = i_code;
                        t[j] = timer_u16(inf_duration.sample(&mut rng));
                        local[n + k] += 1;
                    }
                } else if s[j] == i_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 {
                        s[j] = r_code;
                        t[j] = timer_u16(imm_duration.sample(&mut rng));   // immunity wanes
                        local[2 * n + k] += 1;
                    }
                } else if s[j] == r_code {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 { s[j] = s_code; local[3 * n + k] += 1; }
                }
            }
            local
        })
        .reduce(|| vec![0i64; 4 * n], |mut a, b| { for k in 0..4 * n { a[k] += b[k]; } a });
    list!(
        waned     = block(&tally, 0, n),
        onset     = block(&tally, 1, n),
        recovered = block(&tally, 2, n),
        waned_r   = block(&tally, 3, n)
    )
}

// ── Generic single-transition kernels ────────────────────────────────────────────
//
// The three step kernels above are tuned bundles for the named menagerie. These two are
// the GENERIC building blocks: one timed `from_state -> to_state` transition each, so you
// can compose any model (or one off the menagerie — a vaccinated `V`, a second infectious
// stage, …) from R without writing Rust. Call them downstream-first (move agents out of a
// timed compartment before moving agents in) and apply the returned per-node counts with
// `move_count`. Like the bundles, they branch on the agent's state and decrement the u16
// timer; `step_timer_expire` goes to an UNTIMED destination, `step_timer_expire_set` draws
// the destination's own timer from a Distribution.

/// Generic timed transition `from_state -> to_state` into an UNTIMED destination.
///
/// For each agent in `from_state`, decrements its u16 `timer`; on expiry the agent moves
/// to `to_state` (timer left at 0). Returns the per-node count of transitions. Compose
/// these (downstream-first) to build models beyond the named menagerie; apply the counts
/// with `move_count`. Generalizes the M→S / R→S / I→S or I→R legs.
///
/// @param state,timer,nodeid,count,n_nodes  As in [step_si()].
/// @param from_state Integer state code an agent must occupy to be eligible.
/// @param to_state   Integer (untimed) state code an agent moves to on expiry.
/// @return An integer vector of per-node transition counts (length `n_nodes`).
/// @export
#[extendr]
fn step_timer_expire(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    n_nodes: i32,
    from_state: i32,
    to_state: i32,
) -> Vec<i32> {
    assert!(count >= 0 && n_nodes >= 0, "`count`/`n_nodes` must be non-negative");
    assert!((0..=255).contains(&from_state) && (0..=255).contains(&to_state),
            "`from_state`/`to_state` must be state codes in [0, 255]");
    let count = count as usize;
    let n = n_nodes as usize;
    let (from, to) = (from_state as u8, to_state as u8);

    let chunk = rng::RNG_CHUNK;
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u16_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|((s, t), node)| {
            let mut local = vec![0i64; n];
            for j in 0..s.len() {
                if s[j] == from {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 { s[j] = to; local[node[j] as usize] += 1; }
                }
            }
            local
        })
        .reduce(|| vec![0i64; n], |mut a, b| { for k in 0..n { a[k] += b[k]; } a });
    tally.iter().map(|&x| x as i32).collect()
}

/// Generic timed transition `from_state -> to_state` into a TIMED destination.
///
/// Like [step_timer_expire()] but on expiry the agent's `timer` is reset to a fresh draw
/// from `duration` (the destination state's own clock — e.g. E→I sets the infectious
/// period, I→R sets a waning-immunity period). Returns per-node transition counts.
///
/// @param state,timer,nodeid,count,n_nodes  As in [step_si()].
/// @param from_state Integer state code an agent must occupy to be eligible.
/// @param to_state   Integer state code an agent moves to on expiry.
/// @param duration   A Distribution for the destination state's timer.
/// @return An integer vector of per-node transition counts (length `n_nodes`).
/// @export
#[extendr]
fn step_timer_expire_set(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    n_nodes: i32,
    from_state: i32,
    to_state: i32,
    duration: &Distribution,
) -> Vec<i32> {
    assert!(count >= 0 && n_nodes >= 0, "`count`/`n_nodes` must be non-negative");
    assert!((0..=255).contains(&from_state) && (0..=255).contains(&to_state),
            "`from_state`/`to_state` must be state codes in [0, 255]");
    let count = count as usize;
    let n = n_nodes as usize;
    let (from, to) = (from_state as u8, to_state as u8);

    let base = rng::next_call_base();
    let chunk = rng::RNG_CHUNK;
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u16_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .enumerate()
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|(((ci, s), t), node)| {
            let mut rng = rng::chunk_rng(base, ci);
            let mut local = vec![0i64; n];
            for j in 0..s.len() {
                if s[j] == from {
                    if t[j] > 0 { t[j] -= 1; }
                    if t[j] == 0 {
                        s[j] = to;
                        t[j] = timer_u16(duration.sample(&mut rng));
                        local[node[j] as usize] += 1;
                    }
                }
            }
            local
        })
        .reduce(|| vec![0i64; n], |mut a, b| { for k in 0..n { a[k] += b[k]; } a });
    tally.iter().map(|&x| x as i32).collect()
}

extendr_module! {
    mod steps;
    fn step_si;
    fn step_sir;
    fn step_sirs;
    fn step_timer_expire;
    fn step_timer_expire_set;
}
