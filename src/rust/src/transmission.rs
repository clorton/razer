// ════════════════════════════════════════════════════════════════════════════
// transmission.rs — force of infection and the S→{E,I} transmission kernels.
//
// `calc_foi` computes the per-node force of infection (with spatial coupling) into a
// report Column. The two transmission kernels then move susceptibles into the disease
// chain. Following the project's "kernels mutate agents, return per-node tallies"
// convention, the transmission kernels DO NOT touch the node census: they set each
// infectee's agent `state` (and u16 `timer`) and RETURN the per-node infection count,
// leaving the caller to apply the S↓ / {E,I}↑ census delta and record incidence for
// whichever compartments its model maintains.
//
//   * `transmission`     — S→`to_state` (E or I), drawing the destination's timer from
//                          a Distribution. Used by every model with a timed I.
//   * `transmission_si`  — S→I with I ABSORBING (no timer). Used by the SI model only.
//
// Parallel across agents (Rayon) with a private per-node tally reduced at the end;
// RNG is thread-local (not R-seedable). See sir.rs history / vitals.rs for the idioms.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::column::Column;
use crate::distributions::Distribution;
use crate::epidemic::STATE_S;

// Round a draw to a whole-tick u16 timer, clamped to [1, u16::MAX].
#[inline]
fn timer_u16(x: f64) -> u16 {
    x.round().max(1.0).min(u16::MAX as f64) as u16
}

// Chunk size that splits `count` agents into ~one block per worker thread.
fn chunk_size(count: usize) -> usize {
    let nthreads = rayon::current_num_threads().max(1);
    ((count + nthreads - 1) / nthreads).max(1)
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

// Copy one column (`slot`) of a 2-D Column into an owned Vec<f64>.
fn read_col(col: &Column, slot: usize, n: usize) -> Vec<f64> {
    col.to_f64()[slot * n..(slot + 1) * n].to_vec()
}

/// Compute the per-node force of infection (FOI) for one tick.
///
/// Computes the frequency-dependent, network-redistributed FOI **rate** into column
/// `tick` of `foi`. The local rate per node is
/// `r[k] = beta[k] * seasonality[k] * infected[k] / population[k]`, redistributed as
/// `foi[k] = r[k] * (1 - sum_j W[k, j]) + sum_i r[i] * W[i, k]`. The transmission
/// kernels turn this rate into a per-tick probability `1 - exp(-foi)`.
///
/// Index conventions: `infected` and `population` are census buffers read at the
/// working column `tick + 1` (the denominator is the current population, which may
/// change with vital dynamics); `beta` and `seasonality` are exogenous modifier grids
/// read at the interval column `tick`; the result is written to `foi[tick]`.
///
/// **Ordering and the effective infectious period.** Whether `infected[tick+1]` is the
/// pre- or post-recovery count depends on where you place this call, and that choice
/// sets the realized infectious period:
/// * For **direct S→I** (SIR), call `calc_foi` BEFORE the recovery step, so an agent is
///   still counted on its recovery tick. Because a directly-infected agent is added to
///   `I` *after* this tally (not counted on its entry tick), counting it on its recovery
///   tick instead yields the full period: `R0 = beta * D`, not `beta * (D - 1)`.
/// * For **SEIR-style** entry (agents arrive in `I` via a step run before this tally,
///   e.g. `step_sir`'s E→I), call `calc_foi` AFTER that step: the new infectious are
///   counted on their entry tick and recoveries excluded — also `beta * D`.
///
/// @param infected   Infectious-count census Column (`nodes$I`), `slice_len == n_nodes`.
/// @param population Per-node population census Column (`nodes$N`); the FOI denominator.
/// @param beta       Transmission-coefficient grid (`n_ticks x n_nodes`, from [values_map()]).
/// @param seasonality Seasonal-modifier grid (`n_ticks x n_nodes`, from [values_map()]).
/// @param network    An `n_nodes x n_nodes` numeric coupling matrix.
/// @param foi        A 2-D f64 Column (`(n_ticks-1) x n_nodes`); column `tick` is overwritten.
/// @param tick       0-based tick index: reads `beta[tick]`/`seasonality[tick]` and
///   `infected[tick+1]`/`population[tick+1]`, writes `foi[tick]`.
/// @return `NULL` (invisibly); the result is written into `foi`.
/// @export
#[extendr]
fn calc_foi(
    infected: &Column,
    population: &Column,
    beta: &Column,
    seasonality: &Column,
    network: Robj,
    foi: &mut Column,
    tick: i32,
) {
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let tick = tick as usize;
    let n = foi.slice_len();
    assert!(tick < foi.n_slices(), "`tick` ({tick}) out of range for foi with {} ticks", foi.n_slices());
    for (name, col) in [("infected", &*infected), ("population", &*population),
                        ("beta", &*beta), ("seasonality", &*seasonality)] {
        assert_eq!(col.slice_len(), n, "`{name}` slice length must equal n_nodes ({n})");
    }
    assert!(tick < beta.n_slices(), "`tick` out of range for `beta`");
    assert!(tick < seasonality.n_slices(), "`tick` out of range for `seasonality`");
    assert!(tick + 1 < infected.n_slices(), "`tick`+1 out of range for `infected`");
    assert!(tick + 1 < population.n_slices(), "`tick`+1 out of range for `population`");

    let b   = read_col(beta, tick, n);
    let s   = read_col(seasonality, tick, n);
    let inf = read_col(infected, tick + 1, n);
    let pop = read_col(population, tick + 1, n);
    let w   = read_square_matrix(&network, n);

    let r: Vec<f64> = (0..n)
        .map(|k| if pop[k] > 0.0 { b[k] * s[k] * inf[k] / pop[k] } else { 0.0 })
        .collect();
    let out: Vec<f64> = (0..n)
        .map(|k| {
            let exported: f64 = (0..n).map(|j| w[k + j * n]).sum();
            let imported: f64 = (0..n).map(|i| r[i] * w[i + k * n]).sum();
            r[k] * (1.0 - exported) + imported
        })
        .collect();

    let start = tick * n;
    foi.as_f64_mut()[start..start + n].copy_from_slice(&out);
}

// Per-node infection probability `1 - exp(-foi)` for column `tick` (one exp() per node).
fn infection_probs(foi: &Column, tick: usize, n: usize) -> Vec<f64> {
    foi.as_f64()[tick * n..(tick + 1) * n]
        .iter()
        .map(|&lambda| 1.0 - (-lambda).exp())
        .collect()
}

/// Stochastic transmission S→`to_state`, returning new infections per node.
///
/// Converts column `tick` of `foi` to a per-node probability `1 - exp(-foi)` (once per
/// node), then for each susceptible agent moves it — with that probability — into
/// `to_state` (`E` or `I`), setting its u16 `timer` from `duration` (the incubation or
/// infectious period). The node census is NOT touched: the per-node count of new
/// infections is RETURNED, and the caller applies the `S` ↓ / `to_state` ↑ delta and
/// records incidence as its model requires (it knows whether this is S→E or S→I).
///
/// @param state    Per-agent `u8` state Column (mutated).
/// @param timer    Per-agent `u16` timer Column (mutated; set from `duration`).
/// @param nodeid   Per-agent `u16` 0-based node-id Column.
/// @param count    Number of active agents to process.
/// @param foi      `n_ticks x n_nodes` f64 FOI Column (from [calc_foi()]); column `tick` read.
/// @param tick     0-based tick index.
/// @param to_state State code new infections enter (`laser_states()[["E"]]` or `[["I"]]`).
/// @param duration A Distribution from which the receiving state's timer is drawn.
/// @return An integer vector of new infections per node (length `n_nodes`).
/// @export
#[extendr]
fn transmission(
    state: &mut Column,
    timer: &mut Column,
    nodeid: &Column,
    count: i32,
    foi: &Column,
    tick: i32,
    to_state: i32,
    duration: &Distribution,
) -> Vec<i32> {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    assert!((0..=255).contains(&to_state), "`to_state` must be a state code in [0, 255], got {to_state}");
    let count = count as usize;
    let tick = tick as usize;
    let to = to_state as u8;

    let n = foi.slice_len();
    assert!(tick < foi.n_slices(), "`tick` out of range for `foi`");
    let p = infection_probs(foi, tick, n);

    let susceptible = STATE_S as u8;
    let chunk = chunk_size(count);
    let st = &mut state.as_u8_mut()[..count];
    let tm = &mut timer.as_u16_mut()[..count];
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
        .reduce(|| vec![0i64; n], |mut a, b| { for k in 0..n { a[k] += b[k]; } a });
    new.iter().map(|&x| x as i32).collect()
}

/// Stochastic transmission S→I into an ABSORBING `I` (the SI model), returning new
/// infections per node.
///
/// Like [transmission()] but the agent enters `I` permanently — no `timer` is set (`I`
/// is terminal in SI). Returns the per-node count of new infections; the caller applies
/// the `S` ↓ / `I` ↑ delta.
///
/// @param state    Per-agent `u8` state Column (mutated; S→I).
/// @param nodeid   Per-agent `u16` 0-based node-id Column.
/// @param count    Number of active agents to process.
/// @param foi      `n_ticks x n_nodes` f64 FOI Column; column `tick` read.
/// @param tick     0-based tick index.
/// @return An integer vector of new infections per node (length `n_nodes`).
/// @export
#[extendr]
fn transmission_si(
    state: &mut Column,
    nodeid: &Column,
    count: i32,
    foi: &Column,
    tick: i32,
) -> Vec<i32> {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let count = count as usize;
    let tick = tick as usize;

    let n = foi.slice_len();
    assert!(tick < foi.n_slices(), "`tick` out of range for `foi`");
    let p = infection_probs(foi, tick, n);

    let susceptible = STATE_S as u8;
    let infectious  = crate::epidemic::STATE_I as u8;
    let chunk = chunk_size(count);
    let st = &mut state.as_u8_mut()[..count];
    let nid = &nodeid.as_u16()[..count];
    let new: Vec<i64> = st
        .par_chunks_mut(chunk)
        .zip(nid.par_chunks(chunk))
        .map(|(s, node)| {
            let mut rng = rand::thread_rng();
            let mut local = vec![0i64; n];
            for j in 0..s.len() {
                if s[j] == susceptible {
                    let k = node[j] as usize;
                    if p[k] > 0.0 && rng.gen::<f64>() < p[k] {
                        s[j] = infectious;          // I is absorbing: no timer
                        local[k] += 1;
                    }
                }
            }
            local
        })
        .reduce(|| vec![0i64; n], |mut a, b| { for k in 0..n { a[k] += b[k]; } a });
    new.iter().map(|&x| x as i32).collect()
}

extendr_module! {
    mod transmission;
    fn calc_foi;
    fn transmission;
    fn transmission_si;
}
