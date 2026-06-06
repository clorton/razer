// ════════════════════════════════════════════════════════════════════════════
// mortality.rs — natural (non-disease) mortality by date of death.
//
// Every agent is assigned a date of death `dod` (an absolute tick, drawn at creation
// from a life table — see the Kaplan–Meier estimator). `mortality` is the per-tick
// kernel that retires agents whose scheduled death has arrived: for each living agent
// with `dod <= tick`, it sets the agent's state to D (deceased) and decrements the
// living-compartment census it was counted in, recording the event in the per-node
// deaths flow. It keeps the M/S/E/I/R node census exactly in sync.
//
// Like the other agent-loop kernels it parallelizes across cores (Rayon) with a
// private per-task node buffer reduced at the end (see the project memory: always
// parallelize step/dynamics work across agents).
//
// Orientation for readers coming from C / C++ / C# / Python: see sir.rs / vitals.rs
// for the par_chunks_mut / map-reduce idioms used here. State codes live in epidemic.rs;
// note the u8 `state` Column stores D (= -1 as i32) as 255, so "alive" is `state != 255`.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use crate::column::Column;
use crate::epidemic::{STATE_S, STATE_E, STATE_I, STATE_R, STATE_M, STATE_D};

/// Apply natural mortality for one tick: retire agents whose date of death has arrived.
///
/// For each of the first `count` agents that is still alive (state != D) and whose
/// `dod` (an absolute tick) is `<= tick`, the agent's `state` is set to D and the
/// living compartment it currently occupies (M, S, E, I, or R) is decremented by one at
/// census column `tick + 1` (the working column the caller has already carried
/// forward). The per-node total of deaths this tick is added to the `deaths` flow at
/// column `tick`. Parallelized with private per-thread node buffers summed at the end.
///
/// @param state   Per-agent `u8` state Column (mutated; the deceased become D = 255).
/// @param dod     Per-agent `u32` date-of-death Column (an absolute tick index).
/// @param nodeid  Per-agent `u16` 0-based node-id Column.
/// @param count   Number of active agents to process.
/// @param m_count,s_count,e_count,i_count,r_count  `n_ticks x n_nodes` i32 census
///   Columns kept in sync (mutated at column `tick + 1`).
/// @param deaths  `(n_ticks-1) x n_nodes` i32 flow Column; column `tick` receives the
///   per-node death counts (added).
/// @param tick    0-based tick index.
/// @return `NULL` (invisibly); the Columns are modified in place.
/// @export
#[extendr]
fn mortality(
    state: &mut Column,
    dod: &Column,
    nodeid: &Column,
    count: i32,
    m_count: &mut Column,
    s_count: &mut Column,
    e_count: &mut Column,
    i_count: &mut Column,
    r_count: &mut Column,
    deaths: &mut Column,
    tick: i32,
) {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count  = count as usize;
    let tick_u = tick as u32;            // for the dod <= tick comparison
    let tick   = tick as usize;

    let n = deaths.slice_len();
    for (nm, c) in [("m", &*m_count), ("s", &*s_count), ("e", &*e_count),
                    ("i", &*i_count), ("r", &*r_count)] {
        assert_eq!(c.slice_len(), n, "`{nm}_count` slice length must equal n_nodes");
    }
    assert!(tick < deaths.n_slices(), "`tick` out of range for `deaths`");
    assert!(tick + 1 < m_count.n_slices(), "`tick`+1 out of range for the census buffers");

    // u8 state codes (D is 255 = STATE_D as u8). Slot order in the per-node tally:
    // 0 = M, 1 = S, 2 = E, 3 = I, 4 = R.
    let d_code = STATE_D as u8;
    let m_code = STATE_M as u8;
    let s_code = STATE_S as u8;
    let e_code = STATE_E as u8;
    let i_code = STATE_I as u8;
    let r_code = STATE_R as u8;

    let nthreads = rayon::current_num_threads().max(1);
    let chunk = ((count + nthreads - 1) / nthreads).max(1);
    let st  = &mut state.as_u8_mut()[..count];
    let dd  = &dod.as_u32()[..count];
    let nid = &nodeid.as_u16()[..count];

    // Parallel sweep: each worker keeps a private 5*n tally (M,S,E,I,R decrements per
    // node), marking the newly dead and recording which compartment they left.
    let tally: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(dd.par_chunks(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|((s_chunk, d_chunk), node_chunk)| {
            let mut local = vec![0i64; 5 * n];
            for j in 0..s_chunk.len() {
                let s = s_chunk[j];
                if s != d_code && d_chunk[j] <= tick_u {
                    let k = node_chunk[j] as usize;
                    let slot = if s == m_code { 0 }
                        else if s == s_code { 1 }
                        else if s == e_code { 2 }
                        else if s == i_code { 3 }
                        else if s == r_code { 4 }
                        else { continue };          // unknown alive code: leave untouched
                    local[slot * n + k] += 1;
                    s_chunk[j] = d_code;            // retire the agent
                }
            }
            local
        })
        .reduce(|| vec![0i64; 5 * n], |mut a, b| {
            for k in 0..5 * n { a[k] += b[k]; }
            a
        });

    // Apply the per-compartment census decrements at column tick+1.
    let dst = (tick + 1) * n;
    { let m = m_count.as_i32_mut(); for k in 0..n { m[dst + k] -= tally[k]         as i32; } }
    { let s = s_count.as_i32_mut(); for k in 0..n { s[dst + k] -= tally[n + k]     as i32; } }
    { let e = e_count.as_i32_mut(); for k in 0..n { e[dst + k] -= tally[2 * n + k] as i32; } }
    { let i = i_count.as_i32_mut(); for k in 0..n { i[dst + k] -= tally[3 * n + k] as i32; } }
    { let r = r_count.as_i32_mut(); for k in 0..n { r[dst + k] -= tally[4 * n + k] as i32; } }

    // Record total deaths per node in the flow report at column `tick`.
    let src = tick * n;
    let dthn = deaths.as_i32_mut();
    for k in 0..n {
        let total = tally[k] + tally[n + k] + tally[2 * n + k] + tally[3 * n + k] + tally[4 * n + k];
        dthn[src + k] += total as i32;
    }
}

extendr_module! {
    mod mortality;
    fn mortality;
}
