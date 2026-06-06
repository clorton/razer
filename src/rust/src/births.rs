// ════════════════════════════════════════════════════════════════════════════
// births.rs — crude-birth-rate births, with newborns entering maternal immunity.
//
// Each tick, every living agent gives birth with the per-node daily probability
// `1 - exp(-rate)` (the crude birth rate as a daily hazard, `cbr / 1000 / 365` via
// values_map). Each newborn activates a RESERVED agent slot (past the live `count`,
// sized by calc_capacity) as state M with a `timer` drawn from `maternal_duration`,
// `dob = tick`, and a Kaplan–Meier `dod`. Following the project convention, `births`
// touches no node census: it RETURNS the new active count and the per-node birth count,
// and the caller adds the births to the M census / birth report it maintains.
//
// Two phases (cf. the old step_births_cbr): a PARALLEL per-node birth tally (read-only
// over the living agents), then a SERIAL activation of the new slots (it advances the
// active count and writes new-agent properties). Excess births beyond capacity are
// silently dropped (calc_capacity is meant to make this rare).
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::column::Column;
use crate::distributions::Distribution;
use crate::kmestimator::KaplanMeierEstimator;
use crate::epidemic::{STATE_M, STATE_D};

#[inline]
fn timer_u16(x: f64) -> u16 {
    x.round().max(1.0).min(u16::MAX as f64) as u16
}

/// Apply crude-birth-rate births for one tick; newborns enter maternal immunity (M).
///
/// For each of the first `count` living agents, draws a birth with per-node probability
/// `1 - exp(-birth_rate[node, tick])`. Each birth activates the next reserved slot as a
/// new agent: state M, `timer` from `maternal_duration`, `dob = tick`, and `dod = tick`
/// plus a Kaplan–Meier age at death (from `km`). Returns `list(count, born)`: the new
/// active count (store it back into `people$count`) and the per-node birth count (add it
/// to the `M` census and a birth report). Capped at the allocated capacity.
///
/// @param state    Per-agent `u8` state Column (capacity-sized; newborn slots set to M).
/// @param timer    Per-agent `u16` timer Column (newborn slots set from `maternal_duration`).
/// @param nodeid   Per-agent `u16` node-id Column (newborn slots set to the parent's node).
/// @param dob      Per-agent `i32` date-of-birth Column (newborn slots set to `tick`).
/// @param dod      Per-agent `u32` date-of-death Column (newborn slots set via `km`).
/// @param count    Current active agent count (the first free slot).
/// @param n_nodes  Number of nodes (the length of the returned `born` vector).
/// @param birth_rate `n_ticks x n_nodes` f64 daily-birth-rate grid; column `tick` is read.
/// @param maternal_duration A Distribution for the newborns' maternal-immunity timer.
/// @param km       A KaplanMeierEstimator giving each newborn its age at death.
/// @param tick     0-based tick index.
/// @return `list(count = <new active count>, born = integer[n_nodes])`.
/// @export
#[extendr]
fn births(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &mut Column,
    dob: &mut Column,
    dod: &mut Column,
    count: i32,
    n_nodes: i32,
    birth_rate: &Column,
    maternal_duration: &Distribution,
    km: &KaplanMeierEstimator,
    tick: i32,
) -> List {
    assert!(count >= 0 && n_nodes >= 0, "`count`/`n_nodes` must be non-negative");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count = count as usize;
    let n = n_nodes as usize;
    let tick = tick as usize;
    let capacity = state.len();

    assert_eq!(birth_rate.slice_len(), n, "`birth_rate` slice length must equal n_nodes");
    assert!(tick < birth_rate.n_slices(), "`tick` out of range for `birth_rate`");

    // Per-node daily birth probability (one exp() per node).
    let rate = &birth_rate.as_f64()[tick * n..(tick + 1) * n];
    let p: Vec<f64> = rate.iter().map(|&r| 1.0 - (-r).exp()).collect();
    let d_code = STATE_D as u8;

    // ── Phase 1: parallel per-node birth tally (read-only over the living agents) ─────
    let births_per_node: Vec<i64> = {
        let st  = &state.as_u8()[..count];
        let nid = &nodeid.as_u16()[..count];
        let nthreads = rayon::current_num_threads().max(1);
        let chunk = ((count + nthreads - 1) / nthreads).max(1);
        st.par_chunks(chunk)
            .zip(nid.par_chunks(chunk))
            .map(|(s_chunk, node_chunk)| {
                let mut rng = rand::thread_rng();
                let mut local = vec![0i64; n];
                for j in 0..s_chunk.len() {
                    if s_chunk[j] != d_code {
                        let k = node_chunk[j] as usize;
                        if p[k] > 0.0 && rng.gen::<f64>() < p[k] { local[k] += 1; }
                    }
                }
                local
            })
            .reduce(|| vec![0i64; n], |mut a, b| { for k in 0..n { a[k] += b[k]; } a })
    };

    // ── Phase 2: serial activation of the reserved slots ──────────────────────────────
    let m_code = STATE_M as u8;
    let st   = state.as_u8_mut();
    let tm   = timer.as_u16_mut();
    let nidm = nodeid.as_u16_mut();
    let dobm = dob.as_i32_mut();
    let dodm = dod.as_u32_mut();
    let mut rng = rand::thread_rng();
    let mut idx = count;
    let mut born = vec![0i32; n];
    for k in 0..n {
        let want = births_per_node[k] as usize;
        let made = want.min(capacity - idx);          // cap at capacity (drop the excess)
        for _ in 0..made {
            st[idx]   = m_code;
            tm[idx]   = timer_u16(maternal_duration.sample(&mut rng));
            nidm[idx] = k as u16;
            dobm[idx] = tick as i32;
            dodm[idx] = (tick as i64 + km.sample_newborn_age_at_death(&mut rng)) as u32;
            idx += 1;
        }
        born[k] = made as i32;
        if idx == capacity { break; }
    }

    list!(count = idx as i32, born = born)
}

extendr_module! {
    mod births;
    fn births;
}
