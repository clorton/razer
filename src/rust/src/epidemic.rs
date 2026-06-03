use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::laser_frame::{LaserFrame, PropData};

// ── State constants ───────────────────────────────────────────────────────────
// -1 for D keeps the alive compartments in the non-negative range, so
// `state >= 0` is an unambiguous "is alive?" test without knowing the max
// compartment code.

pub const STATE_S: i32 =  0;
pub const STATE_E: i32 =  1;
pub const STATE_I: i32 =  2;
pub const STATE_R: i32 =  3;
pub const STATE_D: i32 = -1;

/// Named integer vector of epidemic compartment state codes.
///
/// Returns `c(S=0L, E=1L, I=2L, R=3L, D=-1L)`. Use these constants to set
/// and test the `state` property on a people frame.
///
/// @return Named integer vector with elements S, E, I, R, D.
/// @export
#[extendr]
fn laser_states() -> Robj {
    let vals: Vec<i32> = vec![STATE_S, STATE_E, STATE_I, STATE_R, STATE_D];
    let mut robj = vals.iter().collect_robj();
    let names: Vec<String> = ["S", "E", "I", "R", "D"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let names_robj = names.iter().map(|s| s.as_str()).collect_robj();
    robj.set_attrib("names", names_robj).expect("set names");
    robj
}

// ── Raw-pointer helpers ───────────────────────────────────────────────────────
//
// Getting two *mut slices from the same HashMap simultaneously is blocked by
// Rust's aliasing rules: a second get_mut() would be a second exclusive borrow
// while the first is still live.  The solution is to extract a raw *mut i32
// pointer — equivalent to a plain C `int*` — before the Rust borrow is
// released at the end of each block.  The caller is then responsible for
// ensuring the derived slices do not alias and that the backing Vec is not
// structurally modified (no HashMap insert/remove) while the slices are live.
//
// Each function returns as soon as the internal borrow ends, so sequential
// calls on the same frame compile without conflict.

fn int_ptr_mut(frame: &mut LaserFrame, name: &str) -> *mut i32 {
    match &mut frame
        .scalars
        .get_mut(name)
        .unwrap_or_else(|| panic!("frame has no scalar property '{name}'"))
        .data
    {
        PropData::Integer(v) => v.as_mut_ptr(), // equivalent to v.data() in C++
        _ => panic!("property '{name}' must be integer"),
    }
}

fn int_ptr(frame: &LaserFrame, name: &str) -> *const i32 {
    match &frame
        .scalars
        .get(name)
        .unwrap_or_else(|| panic!("frame has no scalar property '{name}'"))
        .data
    {
        PropData::Integer(v) => v.as_ptr(),
        _ => panic!("property '{name}' must be integer"),
    }
}

// ── Transmission kernels ──────────────────────────────────────────────────────

/// Stochastic S→I transmission step (SI / SIR kernel).
///
/// **Parallelism:** the active agent array is split into fixed chunks — one per
/// Rayon worker thread — matching the `prange` pattern. The aggregation phase
/// uses a thread-local fold then reduces by element-wise summation. The FOI
/// phase runs each chunk independently; RNG is thread-local (`rand::thread_rng()`).
///
/// All operations are performed in-place on the backing arrays; no copies are made.
///
/// Node-level I counts are computed from current agent states and written to
/// `nodes$I` (overwriting the previous value).
///
/// **Required people properties:** `state` (integer), `node` (integer, 0-based),
/// `timer` (integer).
/// **Required nodes properties:** `N` (integer, total population per node),
/// `I` (integer, will be overwritten).
///
/// @param people       LaserFrame of agents.
/// @param nodes        LaserFrame of patches/nodes.
/// @param beta         Transmission rate (force of infection per infectious contact per tick).
/// @param inf_duration Infectious period in ticks; written to `timer` on infection.
/// @export
#[extendr]
fn step_transmission_si(
    people: &mut LaserFrame,
    nodes: &mut LaserFrame,
    beta: f64,
    inf_duration: i32,
) {
    let count   = people.count;
    let n_nodes = nodes.count;

    // Extract one raw pointer per property before any slices are created.
    // Each int_ptr_mut() call borrows `people` or `nodes` exclusively for its
    // own duration, then releases it — so sequential calls compile fine.
    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");
    let node_ptr  = int_ptr(people, "node");    // read-only: *const
    let n_pop_ptr = int_ptr(nodes, "N");
    let i_ptr     = int_ptr_mut(nodes, "I");

    // SAFETY: all five pointers come from distinct Vec<i32> allocations.
    // `state` and `timer` are in people.scalars["state"/"timer"]; `node` is in
    // people.scalars["node"]; `n_pop` and `i_out` are in nodes.scalars["N"/"I"].
    // No HashMap entry is inserted or removed while these slices are live, so
    // the Vecs cannot be reallocated or moved.
    // Rayon's par_iter_mut assigns each element to exactly one thread, so
    // concurrent writes to `state` and `timer` do not race.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };
    let node  = unsafe { std::slice::from_raw_parts(node_ptr,  count) };
    let n_pop = unsafe { std::slice::from_raw_parts(n_pop_ptr, n_nodes) };
    let i_out = unsafe { std::slice::from_raw_parts_mut(i_ptr, n_nodes) };

    // ── Parallel aggregation ──────────────────────────────────────────────────
    // fold() is called once per Rayon chunk; each chunk gets its own freshly
    // allocated local Vec<i32> of length n_nodes (the `|| vec![0i32; n_nodes]`
    // closure is the "identity" constructor, analogous to a zero-initialised
    // array in C).  Within a chunk the closure accumulates counts serially.
    // reduce() then merges the per-chunk Vecs pairwise with element-wise
    // addition — equivalent to an OpenMP `reduction(+:i_counts)` clause.
    let i_counts: Vec<i32> = state
        .par_iter()             // parallel read-only iterator over state[0..count]
        .zip(node.par_iter())   // zip with node[0..count]; each thread gets a matching pair of sub-slices
        .fold(
            || vec![0i32; n_nodes],
            |mut local, (&s, &n_idx)| {  // `&s` dereferences the &i32 reference yielded by par_iter
                if s == STATE_I {
                    local[n_idx as usize] += 1;
                }
                local
            },
        )
        .reduce(
            || vec![0i32; n_nodes],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(b.iter()) {
                    *x += y;
                }
                a
            },
        );

    // Write pre-FOI I counts to the nodes frame.
    i_out.copy_from_slice(&i_counts); // equivalent to memcpy

    // ── Parallel stochastic FOI ───────────────────────────────────────────────
    // par_iter_mut() splits state[0..count] into contiguous chunks — one per
    // Rayon worker thread — matching numba's prange pattern.  zip() binds each
    // (s, t) pair with the corresponding n_idx; the compiler sees them as a
    // single unit owned by one thread.
    // rand::thread_rng() is stored in thread-local storage; calling it inside
    // the closure costs only a pointer dereference, with no locking.
    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .zip(node.par_iter())
        .for_each(|((s, t), &n_idx)| {
            if *s == STATE_S {
                let n = n_pop[n_idx as usize] as f64;
                if n > 0.0 {
                    let lambda = beta * i_counts[n_idx as usize] as f64 / n;
                    let p = 1.0 - (-lambda).exp();
                    if rand::thread_rng().gen::<f64>() < p {
                        *s = STATE_I;
                        *t = inf_duration;
                    }
                }
            }
        });
    // Modifications to state[] and timer[] are already in the backing store;
    // no write-back is required.
}

/// Stochastic S→E exposure step (SEIR kernel).
///
/// Same FOI computation and parallelism as `step_transmission_si`, but newly
/// exposed agents move to state E and `timer` is set to `exp_duration`
/// (incubation period). Pair with `step_exposed_ei` to complete E→I.
///
/// **Required people properties:** `state`, `node`, `timer` (all integer).
/// **Required nodes properties:** `N`, `I` (integer; `I` is overwritten).
///
/// @param people        LaserFrame of agents.
/// @param nodes         LaserFrame of patches/nodes.
/// @param beta          Transmission rate.
/// @param exp_duration  Incubation period in ticks; written to `timer` on exposure.
/// @export
#[extendr]
fn step_transmission_se(
    people: &mut LaserFrame,
    nodes: &mut LaserFrame,
    beta: f64,
    exp_duration: i32,
) {
    let count   = people.count;
    let n_nodes = nodes.count;

    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");
    let node_ptr  = int_ptr(people, "node");
    let n_pop_ptr = int_ptr(nodes, "N");
    let i_ptr     = int_ptr_mut(nodes, "I");

    // SAFETY: same invariants as step_transmission_si — five distinct allocations,
    // no HashMap modification while slices are live.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };
    let node  = unsafe { std::slice::from_raw_parts(node_ptr,  count) };
    let n_pop = unsafe { std::slice::from_raw_parts(n_pop_ptr, n_nodes) };
    let i_out = unsafe { std::slice::from_raw_parts_mut(i_ptr, n_nodes) };

    let i_counts: Vec<i32> = state
        .par_iter()
        .zip(node.par_iter())
        .fold(
            || vec![0i32; n_nodes],
            |mut local, (&s, &n_idx)| {
                if s == STATE_I {
                    local[n_idx as usize] += 1;
                }
                local
            },
        )
        .reduce(
            || vec![0i32; n_nodes],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(b.iter()) {
                    *x += y;
                }
                a
            },
        );

    i_out.copy_from_slice(&i_counts);

    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .zip(node.par_iter())
        .for_each(|((s, t), &n_idx)| {
            if *s == STATE_S {
                let n = n_pop[n_idx as usize] as f64;
                if n > 0.0 {
                    let lambda = beta * i_counts[n_idx as usize] as f64 / n;
                    let p = 1.0 - (-lambda).exp();
                    if rand::thread_rng().gen::<f64>() < p {
                        *s = STATE_E;
                        *t = exp_duration;
                    }
                }
            }
        });
}

// ── Timer-based state transitions ─────────────────────────────────────────────
//
// Each of these functions is a pure per-agent operation with no inter-agent
// dependencies, so the Rayon chunk pattern is both correct and optimal.
// All writes are in-place; there is no copy-in / copy-out.

/// Timer-based E→I transition (SEIR kernel).
///
/// For each exposed agent, decrements `timer` by 1. When `timer` reaches 0
/// the agent transitions to state I and `timer` is reset to `inf_duration`.
///
/// @param people        LaserFrame of agents.
/// @param inf_duration  Infectious period in ticks; written to `timer` on E→I.
/// @export
#[extendr]
fn step_exposed_ei(people: &mut LaserFrame, inf_duration: i32) {
    let count     = people.count;
    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");

    // SAFETY: state and timer are distinct Vec allocations; Rayon gives each
    // (s, t) pair to exactly one thread.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };

    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .for_each(|(s, t)| {
            if *s == STATE_E {
                *t -= 1;
                if *t <= 0 {
                    *s = STATE_I;
                    *t = inf_duration;
                }
            }
        });
}

/// Timer-based I→R transition (SIR / SEIR / SEIRS kernel).
///
/// For each infectious agent, decrements `timer` by 1. When `timer` reaches 0
/// the agent transitions to state R and `timer` is reset to `imm_duration`.
///
/// Setting `imm_duration = 0` is correct for SIR and SEIR (where
/// `step_recovered_rs` is never called and the timer value after recovery does
/// not matter). For SEIRS, pass the desired immunity period so that
/// `step_recovered_rs` counts down from the correct starting value.
///
/// @param people        LaserFrame of agents.
/// @param imm_duration  Immunity period in ticks; written to `timer` on I→R.
///   Pass `0L` for SIR / SEIR models where waning is not modelled.
/// @export
#[extendr]
fn step_infectious_ir(people: &mut LaserFrame, imm_duration: i32) {
    let count     = people.count;
    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");

    // SAFETY: same as step_exposed_ei.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };

    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .for_each(|(s, t)| {
            if *s == STATE_I {
                *t -= 1;
                if *t <= 0 {
                    *s = STATE_R;
                    *t = imm_duration; // 0 for SIR/SEIR; imm_duration for SEIRS
                }
            }
        });
}

/// Timer-based I→S transition for SIS models (no immunity).
///
/// For each infectious agent, decrements `timer` by 1. When `timer` reaches 0
/// the agent transitions directly back to state S.
///
/// @param people LaserFrame of agents.
/// @export
#[extendr]
fn step_infectious_is(people: &mut LaserFrame) {
    let count     = people.count;
    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");

    // SAFETY: same as step_exposed_ei.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };

    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .for_each(|(s, t)| {
            if *s == STATE_I {
                *t -= 1;
                if *t <= 0 {
                    *s = STATE_S;
                    *t = 0;
                }
            }
        });
}

/// Timer-based R→S transition for waning-immunity models.
///
/// For each recovered agent, decrements `timer` by 1. When `timer` reaches 0
/// the agent becomes susceptible again (state S).
///
/// @param people LaserFrame of agents.
/// @export
#[extendr]
fn step_recovered_rs(people: &mut LaserFrame) {
    let count     = people.count;
    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");

    // SAFETY: same as step_exposed_ei.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };

    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .for_each(|(s, t)| {
            if *s == STATE_R {
                *t -= 1;
                if *t <= 0 {
                    *s = STATE_S;
                    *t = 0;
                }
            }
        });
}

// ── Demography ────────────────────────────────────────────────────────────────

/// Stochastic mortality step using a crude death rate.
///
/// For each living (non-D) agent, draws Bernoulli(`cdr`) for death.
/// Agents that die have their state set to D (-1) and `timer` set to 0.
/// Dead agents remain in the frame — call `$squash(people$state >= 0L)`
/// periodically to compact.
///
/// Each agent's draw is independent; Rayon assigns a fixed slice to each worker
/// thread. RNG is thread-local (no locking).
///
/// @param people LaserFrame of agents.
/// @param cdr    Crude death rate per agent per tick (probability in \[0, 1\]).
/// @export
#[extendr]
fn step_mortality_cdr(people: &mut LaserFrame, cdr: f64) {
    let count     = people.count;
    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");

    // SAFETY: same as step_exposed_ei.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };

    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .for_each(|(s, t)| {
            if *s != STATE_D && rand::thread_rng().gen::<f64>() < cdr {
                *s = STATE_D;
                *t = 0;
            }
        });
}

/// Stochastic birth step using a crude birth rate.
///
/// **Parallelism:** the Bernoulli draw (does this agent give birth?) is
/// parallelised across all active agents. The subsequent bookkeeping —
/// incrementing `count` and writing new-agent properties — is serial because
/// it modifies the frame's active count.
///
/// For each active (non-D) agent, draws Bernoulli(`cbr`) for a birth event.
/// Each birth creates one new agent that inherits the parent's `node` and
/// starts in state S with `timer = 0`. Other scalar properties of new agents
/// take the default value set at `add_scalar_property` time.
///
/// Excess births beyond `capacity - count` are silently dropped.
///
/// **Required people properties:** `state`, `node`, `timer` (all integer).
///
/// @param people LaserFrame of agents.
/// @param cbr    Crude birth rate per agent per tick (probability in \[0, 1\]).
/// @export
#[extendr]
fn step_births_cbr(people: &mut LaserFrame, cbr: f64) {
    let count    = people.count;
    let capacity = people.capacity;

    // ── Phase 1: parallel Bernoulli draws (read-only) ─────────────────────────
    let state_ptr = int_ptr(people, "state"); // *const — read only
    let node_ptr  = int_ptr(people, "node");

    // SAFETY: const slices; no mutable access to these properties while the
    // slices are live.  Both slices are dropped at the end of this block.
    let (birth_parents, parent_nodes): (Vec<usize>, Vec<i32>) = {
        let state = unsafe { std::slice::from_raw_parts(state_ptr, count) };
        let node  = unsafe { std::slice::from_raw_parts(node_ptr,  count) };

        // filter_map is a combined map+filter: returns Some(idx) to keep, None
        // to discard — equivalent to a loop with `continue` that also transforms.
        let parents: Vec<usize> = state
            .par_iter()
            .enumerate() // yields (index, &element) pairs, like a ranged for with index
            .filter_map(|(idx, &s)| {
                if s != STATE_D && rand::thread_rng().gen::<f64>() < cbr {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect();

        let n = parents.len().min(capacity - count);
        // Collect parent node values before the const slices are dropped.
        let pnodes: Vec<i32> = parents[..n].iter().map(|&idx| node[idx]).collect();
        (parents, pnodes)
        // `state` and `node` slices are dropped here — safe to take *mut next.
    };

    let n_births = parent_nodes.len();
    if n_births == 0 {
        return;
    }

    // ── Phase 2: activate new agents (serial write) ───────────────────────────
    let start_idx = people.count;
    people.count += n_births;

    // Get mutable pointers after the Phase 1 const slices have been dropped.
    let state_ptr = int_ptr_mut(people, "state");
    let node_ptr  = int_ptr_mut(people, "node");
    let timer_ptr = int_ptr_mut(people, "timer");

    // SAFETY: the new-agent range [start_idx, start_idx+n_births) does not
    // overlap the Phase 1 range [0, count).  The backing Vecs were allocated
    // to `capacity` at construction, so this range is within bounds.
    // ptr::add(n) is equivalent to `ptr + n` in C pointer arithmetic.
    unsafe {
        let new_states = std::slice::from_raw_parts_mut(state_ptr.add(start_idx), n_births);
        let new_nodes  = std::slice::from_raw_parts_mut(node_ptr.add(start_idx),  n_births);
        let new_timers = std::slice::from_raw_parts_mut(timer_ptr.add(start_idx), n_births);

        for i in 0..n_births {
            new_states[i] = STATE_S;
            new_nodes[i]  = parent_nodes[i];
            new_timers[i] = 0;
        }
    }

    let _ = birth_parents; // suppress unused-variable warning; len was used above
}

extendr_module! {
    mod epidemic;
    fn laser_states;
    fn step_transmission_si;
    fn step_transmission_se;
    fn step_exposed_ei;
    fn step_infectious_ir;
    fn step_infectious_is;
    fn step_recovered_rs;
    fn step_mortality_cdr;
    fn step_births_cbr;
}
