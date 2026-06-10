# Compartmental (cohort) modeling in Razer: analysis and proposal

*How `laser-cohorts`-style metapopulation/compartmental models map onto
Razer, what is already supported, the precise architectural gaps, and a
proportionate plan to close them. An internal design note (not part of
the published site).*

------------------------------------------------------------------------

## 1. What `laser-cohorts` is

`laser-cohorts` is a **metapopulation / compartmental** model: there are
no agents — only **counts per node per state**, advanced each tick by
binomial draws on difference equations. Its core pieces:

- **`StateArray`** — a NumPy `ndarray` subclass of shape
  `(nticks, n_states, n_nodes)` with a designated `state_axis` (= 1).
  Named access (`states.S`, `states.I`) returns a **zero-copy view**
  along the state axis; `states.S[t]` is the per-node vector at tick
  `t`. Assignment (`states.S[t] -= x`) writes straight into the buffer.
- **`Model`** — holds the `StateArray`, a `network` (`n×n`), parameters,
  and an ordered list of **components**. `run()` carries each state
  forward (`nxt[mask] = cur[mask]`, i.e. copy tick `t → t+1`) then calls
  each component’s `start_step` / `step` / `end_step`.
- **Components** (`Susceptible`, `Infectious`, `InfectiousToRecovered`,
  `TransmissionSI`, `RecoveredToSusceptible`, …) — each owns one
  compartment or transition and mutates the count slices.

The representative dynamics, `TransmissionCommon.step`:

``` python
S = states.S[tick+1]; I = states.I[tick+1]
N = np.maximum(states[tick+1].sum(axis=state_axis-1), 1)   # sum over states -> per-node N
foi[:] = beta[tick] * seasonality[tick] * I / N
transfer = foi[:,None] * network                            # row i scaled by foi[i]
foi += transfer.sum(axis=0)                                 # + inflow
foi -= transfer.sum(axis=1)                                 # - outflow
foi = -np.expm1(-foi)                                       # -> probability
newly_infected = np.random.binomial(S, foi)                 # vectorized draw
S -= newly_infected; sink[tick+1] += newly_infected
```

The operations used: named-state slice views, in-place `+=`/`-=` on
slices, elementwise arithmetic + broadcasting, axis reductions (sum over
states; network row/col sums), `expm1`, and a **vectorized binomial
draw**. Crucially, every array here is **per-node**
(`n_states × n_nodes`), not per-agent.

## 2. How it maps onto Razer

| `laser-cohorts` | Razer today |
|----|----|
| `StateArray` `(nticks, n_states, n_nodes)` | per-state 2-D Columns `allocate_vector(dtype, nticks, n_nodes)`, one per compartment — exactly what `run_model` already allocates (`nodes$S`, `nodes$I`, …) |
| carry-forward `nxt[mask] = cur[mask]` | `carry_forward` / `carry_forward_states` (Rust) |
| `model.network` (`n×n`) | `model$network` (already a 2-D Column) |
| components (`setup`/`start_step`/`step`/`end_step`) | the `init` / `step_enter` / `step_update` / `step_exit` callbacks |
| `states.S[t] -= x` (in-place named slice view) | `nodes$S$col(t)` (read **copy**) + `$set_col(t, …)` (write back) |
| `np.random.binomial`, `expm1`, `sum(axis)`, broadcasting | `rbinom`, `expm1`, `colSums`/`rowSums`, base-R recycling |

## 3. Key finding — it already works today (empirically verified)

A `laser-cohorts`-style spatial SIR (`TransmissionSI` + recovery,
including the network coupling `foi += inflow − outflow`) was
implemented as a 3-node model **entirely on Razer’s existing Columns +
base R** and run for 160 ticks. Result: population conserved per node,
node 0 epidemic (S 9990 → 489), coupling spread it to node 1 (9597
infections from a zero seed), and the uncoupled node 2 stayed at exactly
0.

**Why it just works: compartment arrays are per-*node*, not
per-*agent*.** The rationale for Rust + Columns in Razer is the agent
arrays — millions of elements, copy-on-modify fatal, JIT-less R far too
slow. Compartment counts are `n_states × n_nodes` per tick — a few
thousand elements even for England & Wales (954 nodes). Base R’s
vectorized `+ - * /`, `colSums`, `%*%`, and `rbinom` **are** the
NumPy-operator layer at that scale; copying a 954-vector out of a
Column, computing, and writing it back is microseconds. So the “missing
NumPy operators” are, for compartmental sizes, an **ergonomics**
problem, not a **throughput** problem.

## 4. The actual architectural gaps (ranked, with honest severity)

1.  **Slice ergonomics — “state by node” (real, R-solvable).** `$col(t)`
    returns a *copy*, so you write
    `s <- nodes$S$col(t); …; nodes$S$set_col(t, s - ni)` instead of the
    in-place `states.S[t] -= ni`. The read/compute/write-back round-trip
    is the friction; no Rust needed — a thin R façade closes it.
2.  **A compartmental runner / component framework (missing, R-level).**
    `run_model` is hard-wired to the *agent* kernels; there is no
    count-based runner or component registry. This is the largest
    *functional* gap — and it is pure R to build.
3.  **Unified reproducible RNG (minor).** A pure-R compartmental model
    uses `rbinom` (R’s RNG, `set.seed`), which *is* reproducible and
    single-threaded — not broken, just a *second* RNG stream alongside
    Razer’s `set_seed`. One seeded Rust draw kernel would unify them.
4.  **N-D Columns / axis reductions (not needed now).** `StateArray` is
    3-D with `sum(axis=state)`; Razer uses separate per-state 2-D
    Columns and sums the slices in R. That is fine until you add an
    **age** dimension *and* scale to large cohort grids — only then does
    a true N-D Column + in-Rust reductions earn its keep.

## 5. Proposed solution — layered and proportionate

- **Layer 0 — storage (exists).** Keep per-state 2-D Columns. Add a
  small R `state_array(names, nticks, n_nodes)` S3 façade bundling them
  so you write `st$S`, `st$I`, with helpers `total(st, t)` (per-node N)
  and `slice(st, t)` (an `n_states × n_nodes` matrix for
  sum-over-states). ~40 lines of R.
- **Layer 1 — operators (base R; optional sugar).** Do the per-node
  arithmetic in base R — the right tool at this size. *Optionally* add
  an `Ops.Column` group generic so `colA + colB` / `col * scalar`
  dispatch to in-place Rust elementwise ops, but treat it as ergonomic
  polish, **not** a NumPy-in-Rust: the perf upside at compartment sizes
  is ~nil and the surface area large.
- **Layer 2 — the one worthwhile Rust primitive: seeded vectorized
  draws.** `rbinom_col(n, p) → counts` (and `rpois_col`) drawing
  per-element through Razer’s existing seeded, thread-independent RNG
  (`rng::`). This unifies reproducibility (one `set_seed` controls agent
  *and* compartment models) and is the only piece base R can’t provide
  natively.
- **Layer 3 — a compartmental runner + components.** A
  `run_compartmental()` (or a `mode = "cohort"` on `run_model`) that
  allocates the state Columns, seeds from the scenario, carries forward,
  and runs an ordered list of count components (`Transmission`,
  `Recovery`, `Waning`, `Births`, …) per tick — the direct analog of
  `laser-cohorts`’ component list. **Reuse `calc_foi`** for the network
  coupling (it already computes `foi += inflow − outflow` over the
  network Column), feeding `rbinom_col`.

## 6. Honest recommendation

Do **not** port NumPy into Rust to support compartmental models — that
would be over-engineering driven by a false analogy to the agent case.
Recommended build order:

1.  an R `state_array` façade + a `run_compartmental()` runner with
    count components (pure R, reusing `calc_foi`);
2.  a seeded `rbinom_col` / `rpois_col` Rust kernel for unified
    reproducibility.

That delivers the full `laser-cohorts` menagerie on Razer for a fraction
of the effort. **Defer** N-D Columns and a Column operator-algebra
until/unless age-structured cohorts push the per-tick arrays large
enough that in-Rust ops actually pay off. Suggested first step: a
working `run_compartmental()`, which proves the design end-to-end the
way the agent `run_model` does.
