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
use crate::distributions::Distribution;
use crate::rng;

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
            let mut rng = rng::chunk_rng(base, ci);   // per-chunk seeded RNG
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

/// Import new infectious cases from a schedule, activating reserved agent slots.
///
/// For the given `tick`, scans the schedule (parallel `sched_tick` / `sched_node` /
/// `sched_count` vectors) and, for every entry whose `sched_tick == tick`, activates
/// `sched_count` new agents in node `sched_node`: each takes the next free slot after
/// the active `count`, is set Infectious with a `timer` drawn from `duration`, and
/// has its `nodeid` set. The new agents must fit in the reserved capacity (the
/// `state`/`timer`/`nodeid` Columns are allocated larger than the initial `count`).
/// Per-node import counts are added to the I census at column `tick + 1` and written
/// to the `importations` flow at column `tick`.
///
/// Returns the new active agent count (the caller stores it back into `people$count`).
/// Sequential — it touches only the handful of imported slots, not all agents.
///
/// @param state   Per-agent `u8` state Column (capacity-sized; imported slots set to I).
/// @param timer   Per-agent `u8` timer Column (imported slots set from `duration`).
/// @param nodeid  Per-agent `u16` node-id Column (imported slots set to their node).
/// @param count   Current active agent count (the first free slot).
/// @param i_count `n_ticks x n_nodes` i32 I census Column (mutated at `tick + 1`).
/// @param importations `(n_ticks-1) x n_nodes` i32 flow Column (set at `tick`).
/// @param sched_tick,sched_node,sched_count  Equal-length integer schedule vectors.
/// @param duration A Distribution for the imported cases' infectious timer.
/// @param tick    0-based tick index.
/// @return The new active agent count (`count` plus the number imported this tick).
/// @export
#[extendr]
fn import_infections(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &mut Column,
    count: i32,
    i_count: &mut Column,
    importations: &mut Column,
    sched_tick: Vec<i32>,
    sched_node: Vec<i32>,
    sched_count: Vec<i32>,
    duration: &Distribution,
    tick: i32,
) -> i32 {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    assert_eq!(sched_node.len(), sched_tick.len(), "schedule vectors must have equal length");
    assert_eq!(sched_count.len(), sched_tick.len(), "schedule vectors must have equal length");
    let mut count = count as usize;
    let t = tick as usize;

    let n = importations.slice_len();
    assert_eq!(i_count.slice_len(), n, "`i_count` slice length must equal n_nodes");
    assert!(t < importations.n_slices(), "`tick` out of range for `importations`");
    assert!(t + 1 < i_count.n_slices(), "`tick`+1 out of range for `i_count`");

    // First pass: total this tick + per-node tally, validating nodes and capacity.
    let mut per_node = vec![0i64; n];
    let mut total = 0usize;
    for e in 0..sched_tick.len() {
        if sched_tick[e] == tick {
            let cnt = sched_count[e];
            assert!(cnt >= 0, "schedule counts must be non-negative");
            let node = sched_node[e];
            assert!(node >= 0 && (node as usize) < n, "schedule node {node} out of range [0, {n})");
            per_node[node as usize] += cnt as i64;
            total += cnt as usize;
        }
    }
    let capacity = state.len();
    assert!(
        count + total <= capacity,
        "importing {total} agents would exceed capacity ({capacity}); active count is {count}"
    );

    // Second pass: activate the reserved slots as infectious agents.
    let infectious = STATE_I as u8;
    let st = state.as_u8_mut();
    let tm = timer.as_u16_mut();
    let nid = nodeid.as_u16_mut();
    let mut rng = rng::single_rng();
    for e in 0..sched_tick.len() {
        if sched_tick[e] == tick {
            let node = sched_node[e] as u16;
            for _ in 0..sched_count[e] as usize {
                st[count] = infectious;
                let d = duration.sample(&mut rng);
                tm[count] = d.round().clamp(1.0, 65535.0) as u16;
                nid[count] = node;
                count += 1;
            }
        }
    }

    // Census: add imports to I at tick+1; flow: record imports at tick.
    let dst = (t + 1) * n;
    let ic = i_count.as_i32_mut();
    for k in 0..n { ic[dst + k] += per_node[k] as i32; }
    let imp = importations.as_i32_mut();
    let src = t * n;
    for k in 0..n { imp[src + k] = per_node[k] as i32; }

    count as i32
}

extendr_module! {
    mod vitals;
    fn constant_pop_vitals_sir;
    fn import_infections;
}
