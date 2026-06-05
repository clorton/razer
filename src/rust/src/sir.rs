// ════════════════════════════════════════════════════════════════════════════
// sir.rs — per-tick SIR step kernels over Column buffers.
//
// The node census (S, I, R counts per node, per tick) is maintained INCREMENTALLY:
// each tick carries the previous tick's counts forward and applies only the deltas,
// so the invariant is always  count[t+1] = count[t] ± delta[t]  — no full re-census
// of agents is ever needed. `sir_step` does the carry-forward (for S, I, and R) and
// then the I→R recovery delta; `sir_transmission` applies the S→I infection delta.
//
// Both kernels parallelize the per-agent work across cores (Rayon). The shared
// hazard — counting events per node — is handled with a PRIVATE per-task node
// buffer (no cross-thread writes, no contention), summed (`reduce`) at the end;
// per-node buffers are tiny next to the agent count. (See the project memory:
// always parallelize step/dynamics work across agents.)
//
// Orientation for readers coming from C / C++ / C# / Python:
//   * `slice.par_chunks_mut(k)` (Rayon) splits a `&mut [T]` into consecutive
//     `k`-element chunks handed to worker threads — like an OpenMP parallel-for
//     over blocks. Disjoint chunks ⇒ no data races, so the borrow checker allows
//     concurrent writes. `.zip(...)` pairs the chunks of several slices elementwise.
//   * `.map(|chunk| local).reduce(id, combine)` is map-reduce: each chunk builds a
//     thread-local accumulator (`local`), then `reduce` folds them pairwise into one.
//   * `slice.copy_within(src_range, dest)` is an in-place memmove within one slice
//     (here: copy column `t` onto column `t+1`).
//   * `x.saturating_sub(1)` subtracts but clamps at 0 (no unsigned wraparound).
//   * `dist.sample(&mut rng)` draws from a Distribution; `&Distribution` is `Sync`,
//     so one shared reference is sampled concurrently across worker threads.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::column::Column;
use crate::epidemic::{STATE_S, STATE_I, STATE_R};
use crate::distributions::Distribution;

// Chunk size that splits `count` agents into ~one block per worker thread.
fn chunk_size(count: usize) -> usize {
    let nthreads = rayon::current_num_threads().max(1);
    ((count + nthreads - 1) / nthreads).max(1)
}

/// One SIR recovery step for tick `tick`: recover expired infectious agents.
///
/// For each of the first `count` agents whose `state` is infectious (`I`),
/// decrements its `u8` `timer`; when the timer reaches zero the agent becomes
/// recovered (`R`). The number of recoveries per node is counted in parallel
/// (private per-thread node buffers, summed at the end) and applied as a delta to
/// column `tick + 1` of the census: `I` down, `R` up; the per-node recovery count is
/// also recorded in column `tick` of the `recoveries` flow report.
///
/// The census is updated IN PLACE at column `tick + 1`, which the caller must have
/// already seeded by calling [carry_forward()] on each census counter for this tick
/// (S, I, R) — `sir_step` itself no longer carries the census forward.
///
/// `i_count` / `r_count` are `n_ticks+1 × n_nodes` i32 census Columns; `recoveries`
/// is an `n_ticks × n_nodes` i32 flow Column.
///
/// @param state      Per-agent `u8` state Column (mutated).
/// @param timer      Per-agent `u8` countdown Column (mutated).
/// @param nodeid     Per-agent `u16` 0-based node-id Column.
/// @param count      Number of active agents to process.
/// @param i_count,r_count  `n_ticks+1 × n_nodes` i32 census Columns (mutated at `tick+1`).
/// @param recoveries `n_ticks × n_nodes` i32 flow Column (mutated at `tick`).
/// @param tick       0-based tick index; the recovery delta is applied to column `tick+1`.
/// @return `NULL` (invisibly); the Columns are modified in place.
/// @export
#[extendr]
fn sir_step(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    i_count: &mut Column,
    r_count: &mut Column,
    recoveries: &mut Column,
    tick: i32,
) {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count = count as usize;
    let tick = tick as usize;

    let n = i_count.slice_len();
    assert_eq!(r_count.slice_len(), n, "I/R census must share n_nodes");
    assert_eq!(recoveries.slice_len(), n, "`recoveries` slice length must equal n_nodes");
    assert!(tick + 1 < i_count.n_slices(), "`tick`+1 out of range for the census buffers");
    assert!(tick < recoveries.n_slices(), "`tick` out of range for `recoveries`");

    // The recovery delta lands in column tick+1 (already carry_forward()ed by the
    // caller); the flow lands in column tick.
    let (src, dst) = (tick * n, (tick + 1) * n);

    // Parallel recovery sweep: transition expired infectious agents and tally
    // recoveries per node in a private buffer per worker, then sum the buffers.
    let infectious = STATE_I as u8;
    let recovered = STATE_R as u8;
    let chunk = chunk_size(count);
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u8_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let rec: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|((s, t), node)| {
            let mut local = vec![0i64; n];
            for j in 0..s.len() {
                if s[j] == infectious {
                    t[j] = t[j].saturating_sub(1);
                    if t[j] == 0 {
                        s[j] = recovered;
                        local[node[j] as usize] += 1;
                    }
                }
            }
            local
        })
        .reduce(|| vec![0i64; n], |mut a, b| {
            for k in 0..n { a[k] += b[k]; }
            a
        });

    // Apply the recovery delta to column tick+1 and record the per-node flow.
    let ic = i_count.as_i32_mut();
    for k in 0..n { ic[dst + k] -= rec[k] as i32; }
    let rc = r_count.as_i32_mut();
    for k in 0..n { rc[dst + k] += rec[k] as i32; }
    let flow = recoveries.as_i32_mut();
    for k in 0..n { flow[src + k] = rec[k] as i32; }   // src == tick*n
}

// Read an R numeric/integer matrix into a column-major Vec<f64> of n*n entries.
// `network[i, j]` lives at index `i + j*n` (R matrices are column-major).
fn read_square_matrix(m: &Robj, n: usize) -> Vec<f64> {
    let data: Vec<f64> = if let Some(s) = m.as_real_slice() {
        s.to_vec()
    } else if let Some(s) = m.as_integer_slice() {
        s.iter().map(|&x| x as f64).collect()
    } else {
        panic!("`network` must be a numeric matrix");
    };
    assert_eq!(
        data.len(), n * n,
        "`network` must be an {n} x {n} matrix ({} elements), got {}",
        n * n, data.len()
    );
    data
}

/// Compute the per-node force of infection (FOI) for one tick.
///
/// Reads the infectious count per node from the POST-RECOVERY census column
/// `tick + 1` of `infected` (the working column this tick, after `carry_forward`
/// and `sir_step` have run — so agents recovering this tick are excluded from the
/// force of infection), and with a transmission coefficient `beta` and a spatial
/// coupling `network`, computes the network-redistributed FOI **rate** into column
/// `tick` of `foi`. The local rate is `r[k] = beta[k] * infected[k]`, redistributed
/// as `foi[k] = r[k] * (1 - sum_j W[k, j]) + sum_i r[i] * W[i, k]`.
/// `sir_transmission()` turns this rate into a per-tick probability `1 - exp(-foi)`.
/// Call it after `sir_step` and just before `sir_transmission`.
///
/// `beta` is a single value (every node) or one per node; for frequency-dependent
/// transmission pass `beta = global_beta / N` per node (folding in the `1/N`).
///
/// @param infected An integer/numeric census Column (e.g. `nodes$I`) with
///   `slice_len == n_nodes`; column `tick + 1` (the post-recovery count) is read.
/// @param beta     A numeric vector of length 1 (broadcast) or `n_nodes`.
/// @param network  An `n_nodes x n_nodes` numeric coupling matrix.
/// @param foi      A 2-D f64 Column (`n_ticks x n_nodes`); column `tick` is overwritten.
/// @param tick     0-based tick index: reads `infected[tick+1]`, writes `foi[tick]`.
/// @return `NULL` (invisibly); the result is written into `foi`.
/// @export
#[extendr]
fn calc_foi(infected: &Column, beta: Vec<f64>, network: Robj, foi: &mut Column, tick: i32) {
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let tick = tick as usize;
    let n = foi.slice_len();
    assert!(tick < foi.n_slices(), "`tick` ({tick}) out of range for foi with {} ticks", foi.n_slices());
    assert_eq!(infected.slice_len(), n, "`infected` slice length must equal n_nodes ({n})");
    assert!(tick + 1 < infected.n_slices(), "`tick`+1 ({}) out of range for `infected`", tick + 1);

    // Read the POST-RECOVERY infectious counts: column tick+1 is this tick's working
    // census column (after carry_forward + sir_step), so agents recovering this tick
    // are already excluded from the force of infection.
    let inf_all = infected.to_f64();
    let read = tick + 1;
    let inf = &inf_all[read * n..(read + 1) * n];

    // beta: a single value broadcast to every node, or one value per node.
    let beta_node: Vec<f64> = match beta.len() {
        1 => vec![beta[0]; n],
        len if len == n => beta,
        len => panic!("`beta` must have length 1 or n_nodes ({n}), got {len}"),
    };

    let w = read_square_matrix(&network, n);

    // Local FOI rate, then redistribute through the network (W[i,j] at w[i + j*n]).
    let r: Vec<f64> = (0..n).map(|k| beta_node[k] * inf[k]).collect();
    let out: Vec<f64> = (0..n)
        .map(|k| {
            let exported: f64 = (0..n).map(|j| w[k + j * n]).sum();         // share leaving node k
            let imported: f64 = (0..n).map(|i| r[i] * w[i + k * n]).sum();  // share arriving at node k
            r[k] * (1.0 - exported) + imported
        })
        .collect();

    let start = tick * n;
    foi.as_f64_mut()[start..start + n].copy_from_slice(&out);
}

/// Transmission step for tick `tick`: move susceptibles into a receiving state.
///
/// Converts column `tick` of `foi` (a rate) to a per-node per-tick infection
/// probability `p[k] = 1 - exp(-foi[k])` — computed ONCE per node, not per agent.
/// Then, for each of the first `count` susceptible (`S`) agents, with probability
/// `p[nodeid]` it moves the agent into `to_state` — `I` for an SIR model or `E` for
/// SEIR — drawing the agent's `timer` for that state from `duration` (rounded to
/// whole ticks, clamped to `[1, 255]`) and tallying the event for its node. The
/// timer is whatever clock the receiving state uses next (incubation for `E`,
/// infectious period for `I`); the caller passes the matching `timer` column and
/// `duration`. New infections per node are counted in parallel (private per-thread
/// buffers, summed at the end) and applied as a delta to column `tick + 1` of the
/// census (`S` down, the receiving state's `to_count` up; the caller has already
/// carried the census forward), and recorded in column `tick` of `incidence`.
///
/// @param state      Per-agent `u8` state Column (mutated).
/// @param timer      Per-agent `u8` countdown Column for the receiving state (mutated).
/// @param nodeid     Per-agent `u16` 0-based node-id Column.
/// @param count      Number of active agents to process.
/// @param foi        `n_ticks x n_nodes` f64 FOI Column (from [calc_foi()]); column
///   `tick` is read.
/// @param s_count    `n_ticks+1 × n_nodes` i32 census Column for `S` (mutated at `tick+1`).
/// @param to_count   `n_ticks+1 × n_nodes` i32 census Column for the receiving state
///   `to_state` (mutated at `tick+1`).
/// @param incidence  `n_ticks x n_nodes` i32 flow Column (mutated at `tick`).
/// @param tick       0-based tick index.
/// @param to_state   State code new infections enter (e.g. `laser_states()[["I"]]`
///   for SIR or `[["E"]]` for SEIR).
/// @param duration   A Distribution from which the receiving state's timer is drawn.
/// @return `NULL` (invisibly); the Columns are modified in place.
/// @export
#[extendr]
fn transmission(
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

    // Per-NODE infection probability (one exp() per node, reused for every agent in
    // that node) — far cheaper than 1 - exp(-foi) per susceptible agent.
    let foi_col = &foi.as_f64()[tick * n..(tick + 1) * n];
    let p: Vec<f64> = foi_col.iter().map(|&lambda| 1.0 - (-lambda).exp()).collect();

    // Parallel transmission sweep: move susceptibles into `to_state` and tally the
    // events per node in a private buffer per worker, then sum the buffers.
    let susceptible = STATE_S as u8;
    let chunk = chunk_size(count);
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u8_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let new: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(tm.par_chunks_mut(chunk))
        .zip(nid.par_chunks(chunk))
        .map(|((s, t), node)| {
            let mut rng = rand::thread_rng();   // per-worker thread-local RNG
            let mut local = vec![0i64; n];
            for j in 0..s.len() {
                if s[j] == susceptible {
                    let k = node[j] as usize;
                    if p[k] > 0.0 && rng.gen::<f64>() < p[k] {
                        s[j] = to;
                        let d = duration.sample(&mut rng);
                        t[j] = d.round().clamp(1.0, 255.0) as u8;
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

    // Apply the delta to column tick+1 (S down, receiving state up) and record the flow.
    let dst = (tick + 1) * n;
    let sc = s_count.as_i32_mut();
    for k in 0..n { sc[dst + k] -= new[k] as i32; }
    let tc = to_count.as_i32_mut();
    for k in 0..n { tc[dst + k] += new[k] as i32; }
    let flow = incidence.as_i32_mut();
    for k in 0..n { flow[tick * n + k] = new[k] as i32; }
}

extendr_module! {
    mod sir;
    fn sir_step;
    fn calc_foi;
    fn transmission;
}
