// ════════════════════════════════════════════════════════════════════════════
// vitals.rs — vital-dynamics kernels over Column buffers.
//
// `constant_pop_vitals_sir` models births and deaths at a constant population for an
// SIR model: each agent dies with a per-tick probability derived from a crude death
// rate, and a death is immediately replaced by a newborn susceptible in the same
// slot (state reset to S, timer reset to 0). Every such event is recorded as BOTH a
// birth and a death (they balance under constant population), and the S/I/R node
// census is kept exactly in sync: an agent dying out of I or R moves that count down
// and S up (a death out of S nets to zero). This resupplies susceptibles and enables
// endemic dynamics. The kernel is SIR-specific (it knows the S, I, R compartments);
// an SEIR variant would also track E.
//
// Like the other agent-loop kernels it parallelizes across cores (Rayon) with a
// private per-task node buffer reduced at the end (see the project memory: always
// parallelize step/dynamics work across agents).
//
// Orientation for readers coming from C / C++ / C# / Python: see sir.rs for the
// par_chunks_mut / map-reduce / thread-local-RNG idioms used here.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::column::Column;
use crate::epidemic::{STATE_S, STATE_I, STATE_R};

/// Apply constant-population SIR vital dynamics for one tick.
///
/// `rate` is a per-node daily death-HAZARD-rate grid (a values map; the caller
/// converts a crude death rate of annual deaths per 1000 people to a daily rate,
/// e.g. `cdr / 1000 / 365`). For each of the first `count` agents, this converts the
/// node's rate to a probability `1 - exp(-rate)` (once per node, like transmission)
/// and, with that probability, "unalives" the agent and replaces it with a newborn:
/// the agent's `state` is reset to Susceptible and its `timer` to 0.
///
/// The S/I/R node census is updated IN PLACE at column `tick + 1` (the working
/// column the caller has already carried forward): a death out of I decrements `I`
/// and increments `S`; a death out of R decrements `R` and increments `S`; a death
/// out of S nets to zero. Every event (from any compartment) is counted per node and
/// written to BOTH the `births` and `deaths` flow reports for `tick` (equal under
/// constant population). Agents are assumed to be in S, I, or R (it is the SIR
/// variant). Parallelized with private per-thread node buffers summed at the end.
///
/// @param state   Per-agent `u8` state Column (mutated; deaths reset to Susceptible).
/// @param timer   Per-agent `u8` countdown Column (mutated; deaths reset to 0).
/// @param nodeid  Per-agent `u16` 0-based node-id Column.
/// @param count   Number of active agents to process.
/// @param rate    Per-node daily death-hazard-rate grid (`n_ticks x n_nodes`, from
///   [values_map()]); column `tick` is read.
/// @param s_count,i_count,r_count  `n_ticks x n_nodes` i32 census Columns kept in
///   sync (mutated at column `tick + 1`).
/// @param births,deaths  `(n_ticks-1) x n_nodes` i32 flow Columns; column `tick`
///   receives the per-node event counts (equal; mutated).
/// @param tick    0-based tick index.
/// @return `NULL` (invisibly); the Columns are modified in place.
/// @export
#[extendr]
fn constant_pop_vitals_sir(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    rate: &Column,
    s_count: &mut Column,
    i_count: &mut Column,
    r_count: &mut Column,
    births: &mut Column,
    deaths: &mut Column,
    tick: i32,
) {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count = count as usize;
    let tick = tick as usize;

    let n = births.slice_len();
    assert_eq!(deaths.slice_len(), n, "`births`/`deaths` must share n_nodes");
    assert_eq!(rate.slice_len(), n, "`rate` slice length must equal n_nodes");
    assert_eq!(s_count.slice_len(), n, "S/I/R census must share n_nodes");
    assert_eq!(i_count.slice_len(), n, "S/I/R census must share n_nodes");
    assert_eq!(r_count.slice_len(), n, "S/I/R census must share n_nodes");
    assert!(tick < rate.n_slices(), "`tick` out of range for `rate`");
    assert!(tick < births.n_slices(), "`tick` out of range for `births`");
    assert!(tick < deaths.n_slices(), "`tick` out of range for `deaths`");
    assert!(tick + 1 < i_count.n_slices(), "`tick`+1 out of range for the census buffers");

    // Per-node death probability for this tick (one exp() per node).
    let rate_all = rate.to_f64();
    let p: Vec<f64> = rate_all[tick * n..(tick + 1) * n]
        .iter()
        .map(|&r| 1.0 - (-r).exp())
        .collect();

    // Parallel sweep. Each worker keeps a private per-node buffer of three tallies
    // (laid out as 3 contiguous n-length blocks): total events [0..n], deaths out of
    // I [n..2n], deaths out of R [2n..3n]. A death out of S contributes to the total
    // only (S -> S nets to zero in the census).
    let s_code = STATE_S as u8;
    let i_code = STATE_I as u8;
    let r_code = STATE_R as u8;
    let nthreads = rayon::current_num_threads().max(1);
    let chunk = ((count + nthreads - 1) / nthreads).max(1);
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u8_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|((s, t), node)| {
            let mut rng = rand::thread_rng();   // per-worker thread-local RNG
            let mut local = vec![0i64; 3 * n];
            for j in 0..s.len() {
                let k = node[j] as usize;
                if p[k] > 0.0 && rng.gen::<f64>() < p[k] {
                    let old = s[j];
                    local[k] += 1;                                    // total events
                    if old == i_code { local[n + k] += 1; }           // death out of I
                    else if old == r_code { local[2 * n + k] += 1; }  // death out of R
                    debug_assert!(old == s_code || old == i_code || old == r_code,
                                  "constant_pop_vitals_sir: agent in a non-SIR state");
                    s[j] = s_code;   // death replaced by a newborn susceptible
                    t[j] = 0;        // reset the agent's timer
                }
            }
            local
        })
        .reduce(|| vec![0i64; 3 * n], |mut a, b| {
            for k in 0..3 * n { a[k] += b[k]; }
            a
        });

    // Apply the per-compartment census delta at column tick+1: S up, I/R down.
    let dst = (tick + 1) * n;
    let ic = i_count.as_i32_mut();
    for k in 0..n { ic[dst + k] -= tally[n + k] as i32; }            // - deaths out of I
    let rc = r_count.as_i32_mut();
    for k in 0..n { rc[dst + k] -= tally[2 * n + k] as i32; }        // - deaths out of R
    let sc = s_count.as_i32_mut();
    for k in 0..n { sc[dst + k] += (tally[n + k] + tally[2 * n + k]) as i32; } // + reborn from I/R

    // Record total events per node as both births and deaths for this tick.
    let start = tick * n;
    let b = births.as_i32_mut();
    for k in 0..n { b[start + k] = tally[k] as i32; }
    let d = deaths.as_i32_mut();
    for k in 0..n { d[start + k] = tally[k] as i32; }
}

extendr_module! {
    mod vitals;
    fn constant_pop_vitals_sir;
}
