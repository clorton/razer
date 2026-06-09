// ════════════════════════════════════════════════════════════════════════════
// mortality.rs — natural (non-disease) mortality by date of death.
//
// Every agent is assigned a date of death `dod` (an absolute tick) at creation. Each
// tick, `mortality` retires the agents whose scheduled death has arrived: it sets their
// state to D (deceased) and RETURNS the per-node count of deaths broken down by the
// state each agent left, so the caller can decrement whichever states its
// model maintains (and total them into a deaths report). It touches no node census.
//
// Parallel across agents (Rayon) with a private per-node tally reduced at the end.
// State codes live in epidemic.rs; the u8 `state` Column stores D (= -1) as 255, so
// "alive" is `state != 255`.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use crate::column::Column;
use crate::epidemic::{STATE_S, STATE_E, STATE_I, STATE_R, STATE_M, STATE_D};

/// Apply natural mortality for one tick, returning deaths per node by state.
///
/// For each of the first `count` living agents (state != D) whose `dod` (an absolute
/// tick) is `<= tick`, sets the agent's `state` to D and tallies the death against the
/// state it occupied. Returns `list(m, s, e, i, r)` of per-node death counts; the
/// caller decrements those census states (and records the total deaths flow).
///
/// @param state   Per-agent `u8` state Column (mutated; the deceased become D = 255).
/// @param dod     Per-agent `u32` date-of-death Column (an absolute tick index).
/// @param nodeid  Per-agent `u16` 0-based node-id Column.
/// @param count   Number of active agents to process.
/// @param n_nodes Number of nodes (the length of each returned vector).
/// @param tick    0-based tick index; agents with `dod <= tick` die.
/// @return `list(m, s, e, i, r)` of `integer[n_nodes]` death counts by source state.
/// @export
#[extendr]
fn mortality(
    state: &mut Column,
    dod: &Column,
    nodeid: &Column,
    count: i32,
    n_nodes: i32,
    tick: i32,
) -> List {
    assert!(count >= 0 && n_nodes >= 0, "`count`/`n_nodes` must be non-negative");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count  = count as usize;
    let n      = n_nodes as usize;
    let tick_u = tick as u32;

    // Slot order in the per-node tally: 0 = M, 1 = S, 2 = E, 3 = I, 4 = R.
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
                        else { continue };
                    local[slot * n + k] += 1;
                    s_chunk[j] = d_code;
                }
            }
            local
        })
        .reduce(|| vec![0i64; 5 * n], |mut a, b| { for k in 0..5 * n { a[k] += b[k]; } a });

    let blk = |slot: usize| -> Vec<i32> {
        tally[slot * n..(slot + 1) * n].iter().map(|&x| x as i32).collect()
    };
    list!(m = blk(0), s = blk(1), e = blk(2), i = blk(3), r = blk(4))
}

extendr_module! {
    mod mortality;
    fn mortality;
}
