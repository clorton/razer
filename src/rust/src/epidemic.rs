use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::laser_frame::{LaserFrame, PropData};
use crate::distributions::Distribution;

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

// ── Duration helper ─────────────────────────────────────────────────────────────
//
// Distributions return floating-point draws; a state timer is a whole number of
// ticks. Each step rounds its draw to the nearest tick and clamps to a minimum of
// 1, so a freshly entered timed state always lasts at least one tick (a timer of 0
// would be decremented to a transition on the very next step). Rounding (rather
// than truncation) keeps the realized mean duration unbiased relative to the
// requested distribution mean.
#[inline]
fn duration_ticks(x: f64) -> i32 {
    x.round().max(1.0) as i32
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
/// @param inf_dist     A `Distribution` (e.g. `dist_constant()` or `dist_normal()`) giving the
///   infectious period in ticks; sampled per newly infected agent and written to
///   `timer` on S→I (rounded to whole ticks, clamped to a minimum of 1).
/// @export
#[extendr]
fn step_transmission_si(
    people: &mut LaserFrame,
    nodes: &mut LaserFrame,
    beta: f64,
    inf_dist: &Distribution,
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
    // for_each_init() builds one thread_rng per worker (Pattern B: the kernel
    // owns the RNG); the same rng drives both the Bernoulli infection draw and
    // the infectious-period draw from the shared `inf_dist`.
    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .zip(node.par_iter())
        .for_each_init(rand::thread_rng, |rng, ((s, t), &n_idx)| {
            if *s == STATE_S {
                let n = n_pop[n_idx as usize] as f64;
                if n > 0.0 {
                    let lambda = beta * i_counts[n_idx as usize] as f64 / n;
                    let p = 1.0 - (-lambda).exp();
                    if rng.gen::<f64>() < p {
                        *s = STATE_I;
                        *t = duration_ticks(inf_dist.sample(rng));
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
/// exposed agents move to state E and `timer` is set to a draw from `exp_dist`
/// (incubation period). Pair with `step_exposed_ei` to complete E→I.
///
/// **Required people properties:** `state`, `node`, `timer` (all integer).
/// **Required nodes properties:** `N`, `I` (integer; `I` is overwritten).
///
/// @param people        LaserFrame of agents.
/// @param nodes         LaserFrame of patches/nodes.
/// @param beta          Transmission rate.
/// @param exp_dist      A `Distribution` (e.g. `dist_constant()` or `dist_normal()`) giving the
///   incubation period in ticks; sampled per newly exposed agent and written to
///   `timer` on S→E (rounded to whole ticks, clamped to a minimum of 1).
/// @export
#[extendr]
fn step_transmission_se(
    people: &mut LaserFrame,
    nodes: &mut LaserFrame,
    beta: f64,
    exp_dist: &Distribution,
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
        .for_each_init(rand::thread_rng, |rng, ((s, t), &n_idx)| {
            if *s == STATE_S {
                let n = n_pop[n_idx as usize] as f64;
                if n > 0.0 {
                    let lambda = beta * i_counts[n_idx as usize] as f64 / n;
                    let p = 1.0 - (-lambda).exp();
                    if rng.gen::<f64>() < p {
                        *s = STATE_E;
                        *t = duration_ticks(exp_dist.sample(rng));
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
//
// Two generalized kernels (mirroring laser-generic's `nb_timer_update` and
// `nb_timer_update_timer_set`) capture the full pattern of timer-driven state
// change, parameterized by the `from`/`to` state codes:
//
//   * `step_timer_expire`     — transition into an *absorbing* (untimed) state;
//                               the timer is left at 0 on arrival.
//   * `step_timer_expire_set` — transition into a state that has *its own*
//                               duration; a fresh per-agent timer is drawn from
//                               a `Distribution` on arrival.
//
// The four named kernels below (`step_exposed_ei`, `step_infectious_ir`,
// `step_infectious_is`, `step_recovered_rs`) are thin, fixed-state wrappers over
// these two helpers, so the model-building vocabulary stays readable while the
// core loop lives in exactly one place.

// Private core loop for the absorbing-state transition. Decrements the timer of
// every agent in `from_state`; on expiry the agent moves to `to_state` with its
// timer reset to 0 (the destination carries no duration of its own).
fn timer_expire_impl(people: &mut LaserFrame, from_state: i32, to_state: i32) {
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
            if *s == from_state {
                *t -= 1;
                if *t <= 0 {
                    *s = to_state;
                    *t = 0;
                }
            }
        });
}

// Private core loop for the timed-destination transition. Like `timer_expire_impl`
// but, on expiry, resets the timer to a fresh draw from `dist` (rounded to whole
// ticks, clamped to a minimum of 1) — the destination state's own duration.
fn timer_expire_set_impl(
    people: &mut LaserFrame,
    from_state: i32,
    to_state: i32,
    dist: &Distribution,
) {
    let count     = people.count;
    let state_ptr = int_ptr_mut(people, "state");
    let timer_ptr = int_ptr_mut(people, "timer");

    // SAFETY: state and timer are distinct Vec allocations; Rayon gives each
    // (s, t) pair to exactly one thread.
    let state = unsafe { std::slice::from_raw_parts_mut(state_ptr, count) };
    let timer = unsafe { std::slice::from_raw_parts_mut(timer_ptr, count) };

    // for_each_init constructs one thread_rng per Rayon worker and threads it
    // through that worker's chunk, so the thread-local lookup is paid once per
    // chunk rather than once per agent.
    state
        .par_iter_mut()
        .zip(timer.par_iter_mut())
        .for_each_init(rand::thread_rng, |rng, (s, t)| {
            if *s == from_state {
                *t -= 1;
                if *t <= 0 {
                    *s = to_state;
                    *t = duration_ticks(dist.sample(rng));
                }
            }
        });
}

/// Generalized timer-expiry transition into an *absorbing* (untimed) state.
///
/// For each agent currently in `from_state`, decrements `timer` by 1; when the
/// timer reaches 0 the agent moves to `to_state` and its timer is left at 0. The
/// destination carries no duration of its own — this is the "transition to an
/// absorbing state" generalization (laser-generic's `nb_timer_update`), e.g.
/// I→S in SIS, I→R in SIR, or R→S waning.
///
/// This is the engine behind [step_infectious_is()] (I→S) and
/// [step_recovered_rs()] (R→S):
/// `step_timer_expire(people, laser_states()[["I"]], laser_states()[["S"]])`
/// is exactly `step_infectious_is(people)`.
///
/// **Required people properties:** `state`, `timer` (both integer).
///
/// @param people     LaserFrame of agents.
/// @param from_state Integer state code an agent must currently occupy to be eligible.
/// @param to_state   Integer state code an agent moves to when its timer expires.
/// @export
#[extendr]
fn step_timer_expire(people: &mut LaserFrame, from_state: i32, to_state: i32) {
    timer_expire_impl(people, from_state, to_state);
}

/// Generalized timer-expiry transition into a state that has *its own* duration.
///
/// For each agent currently in `from_state`, decrements `timer` by 1; when the
/// timer reaches 0 the agent moves to `to_state` and `timer` is reset to a fresh
/// per-agent draw from `duration_dist` (rounded to whole ticks, clamped to a
/// minimum of 1). This is the "transition to a state with its own duration timer"
/// generalization (laser-generic's `nb_timer_update_timer_set`), e.g. E→I in
/// SEIR or I→R with waning in SEIRS.
///
/// This is the engine behind [step_exposed_ei()] (E→I) and
/// [step_infectious_ir()] (I→R):
/// `step_timer_expire_set(people, laser_states()[["E"]], laser_states()[["I"]], inf_dist)`
/// is exactly `step_exposed_ei(people, inf_dist)`.
///
/// **RNG:** thread-local (Pattern B) — each Rayon worker draws from its own
/// `thread_rng`; the single `duration_dist` handle is shared across threads by
/// reference.
///
/// **Required people properties:** `state`, `timer` (both integer).
///
/// @param people        LaserFrame of agents.
/// @param from_state    Integer state code an agent must currently occupy to be eligible.
/// @param to_state      Integer state code an agent moves to when its timer expires.
/// @param duration_dist A `Distribution` (e.g. `dist_constant()` or `dist_normal()`)
///   giving the destination state's duration in ticks; sampled per transitioning
///   agent and written to `timer`.
/// @export
#[extendr]
fn step_timer_expire_set(
    people: &mut LaserFrame,
    from_state: i32,
    to_state: i32,
    duration_dist: &Distribution,
) {
    timer_expire_set_impl(people, from_state, to_state, duration_dist);
}

/// Timer-based E→I transition (SEIR kernel).
///
/// For each exposed agent, decrements `timer` by 1. When `timer` reaches 0 the
/// agent transitions to state I and `timer` is set to a fresh draw from
/// `inf_dist`, the infectious-period distribution. Pass `dist_constant(d)` for a
/// fixed period of `d` ticks, or e.g. `dist_normal(mean, variance)` for a stochastic
/// per-agent period. Draws are rounded to the nearest tick and clamped to a
/// minimum of 1.
///
/// A fixed-state shorthand for
/// [step_timer_expire_set()]`(people, E, I, inf_dist)`.
///
/// **RNG:** thread-local — each Rayon worker draws from its own `thread_rng`
/// (Pattern B: the kernel owns the RNG and passes it into the sampler). The
/// single `inf_dist` handle is shared across threads by reference.
///
/// @param people    LaserFrame of agents.
/// @param inf_dist  A `Distribution` (e.g. from `dist_constant()` or `dist_normal()`)
///   giving the infectious period in ticks; sampled and written to `timer` on E→I.
/// @export
#[extendr]
fn step_exposed_ei(people: &mut LaserFrame, inf_dist: &Distribution) {
    timer_expire_set_impl(people, STATE_E, STATE_I, inf_dist);
}

/// Timer-based I→R transition (SIR / SEIR / SEIRS kernel).
///
/// For each infectious agent, decrements `timer` by 1. When `timer` reaches 0
/// the agent transitions to state R and `timer` is set to a draw from `imm_dist`,
/// the immunity (waning) period.
///
/// A fixed-state shorthand for
/// [step_timer_expire_set()]`(people, I, R, imm_dist)`.
///
/// For SEIRS, pass the desired immunity distribution so that `step_recovered_rs`
/// counts down R→S from a fresh per-agent draw. For SIR and SEIR (no waning,
/// `step_recovered_rs` never called) the R timer is never read, so any
/// distribution works — `dist_constant(0)` is the conventional "don't care".
///
/// **RNG:** thread-local — each Rayon worker draws from its own `thread_rng`
/// (Pattern B). The single `imm_dist` handle is shared across threads by reference.
///
/// @param people    LaserFrame of agents.
/// @param imm_dist  A `Distribution` (e.g. `dist_constant()` or `dist_normal()`) giving the
///   immunity period in ticks; sampled and written to `timer` on I→R (rounded to
///   whole ticks, clamped to a minimum of 1). Use `dist_constant(0)` for SIR / SEIR.
/// @export
#[extendr]
fn step_infectious_ir(people: &mut LaserFrame, imm_dist: &Distribution) {
    timer_expire_set_impl(people, STATE_I, STATE_R, imm_dist);
}

/// Timer-based I→S transition for SIS models (no immunity).
///
/// For each infectious agent, decrements `timer` by 1. When `timer` reaches 0
/// the agent transitions directly back to state S.
///
/// A fixed-state shorthand for [step_timer_expire()]`(people, I, S)`.
///
/// @param people LaserFrame of agents.
/// @export
#[extendr]
fn step_infectious_is(people: &mut LaserFrame) {
    timer_expire_impl(people, STATE_I, STATE_S);
}

/// Timer-based R→S transition for waning-immunity models.
///
/// For each recovered agent, decrements `timer` by 1. When `timer` reaches 0
/// the agent becomes susceptible again (state S).
///
/// A fixed-state shorthand for [step_timer_expire()]`(people, R, S)`.
///
/// @param people LaserFrame of agents.
/// @export
#[extendr]
fn step_recovered_rs(people: &mut LaserFrame) {
    timer_expire_impl(people, STATE_R, STATE_S);
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
    fn step_timer_expire;
    fn step_timer_expire_set;
    fn step_exposed_ei;
    fn step_infectious_ir;
    fn step_infectious_is;
    fn step_recovered_rs;
    fn step_mortality_cdr;
    fn step_births_cbr;
}
