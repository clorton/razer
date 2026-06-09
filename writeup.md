# Bringing the LASER architecture to R: an assessment of `razer`

*An honest engineering retrospective comparing the Python **LASER** stack
(`laser-core` + `laser-generic`) with the R/Rust port, **`razer`** — what carried over,
what didn't, what's genuinely better, and whether an R user should reach for `razer` or
just drive LASER through `reticulate`.*

---

## 0. Scope and sources

This document compares two living codebases:

- **LASER** — the Institute for Disease Modeling's Python toolkit for fast, composable,
  spatial agent-based disease models. Two layers matter here: **`laser-core`** (the
  performance substrate: the `LaserFrame` struct-of-arrays container, `PropertySet`,
  distance/migration models, demographic samplers, the Kaplan–Meier estimator, RNG, and
  Numba-compiled kernels) and **`laser-generic`** (a component library that assembles
  SI/SIR/SEIR-family models out of reusable processes).
- **`razer`** — the R package in this repository, an extendr/Rust port of that
  architecture.

The `razer` side of every claim is grounded in this repository's source. The LASER side is
grounded in its **documented architecture** (the `software-overview` and tutorial docs) and
its established public design; where an exact Python API signature would matter and I
cannot verify it against source in this tree, I describe the behavior rather than invent a
signature. Treat the LASER specifics as architecturally accurate but not line-cited.

---

## 1. What LASER is, and why it is fast

LASER's core insight is a **struct-of-arrays (SoA)** memory layout. Instead of one Python
object per agent, every agent *property* is a flat array and an active population is a
contiguous prefix of that array. This is the `LaserFrame`: a capacity-sized set of parallel
NumPy arrays (`state`, `nodeid`, timers, dates of birth/death, …) plus a live `count`,
with helpers to add scalar/vector properties, grow the population into reserved capacity,
and `squash` the dead out so slots can be reused. A second `LaserFrame` typically holds the
per-node / per-patch report time series.

Performance comes from a **two-tool stack**:

1. **NumPy** supplies the typed, contiguous arrays and all the vectorized bulk operations
   (reductions, masked assignment, slicing) — cache-friendly and implemented in C.
2. **Numba** JIT-compiles the *hot per-agent loops* — transmission, timer decrements, births,
   mortality — to native code, with `@njit(parallel=True)` and `prange` giving thread-level
   parallelism over the agent array.

`laser-generic` then layers a **component model** on top: a `Model` owns the `people` and
`patches` `LaserFrame`s and a `PropertySet` of parameters, and runs an ordered list of
**components** (e.g. `Births`, `Mortality`, `Susceptibility`, `Transmission`, the
exposed/infectious/recovered progressions) once per tick. Each component is a class with
initialization and per-tick `step` behavior; a model is *composed* by choosing which
components to register. New science = a new component class, no core changes.

The whole thing is **Python-native**: the arrays are NumPy arrays the user can slice, plot,
and analyze directly; parameters are attribute-accessible; components are introspectable
objects.

---

## 2. The R problem: no NumPy, no Numba

R cannot replicate that stack directly, for two structural reasons:

- **No NumPy equivalent.** Base R vectors are limited to `logical`, `integer` (32-bit
  signed), `double`, and `complex` — there are **no unsigned or narrow integer types**
  (`u8`/`u16`/`u32`), which matter enormously at national/global agent counts (a `u8` state
  \+ `u16` node id + `u16` timer is 5 bytes/agent; the same in R doubles is 24). Worse, R has
  **copy-on-modify** semantics: mutating an array in a function generally copies it. There is
  no in-place, mutable, typed buffer that hot loops can hammer without allocation.
- **No Numba equivalent.** R has no production JIT for tight numerical loops. The byte-code
  compiler does not approach native speed, and a per-agent `for` loop over millions of
  agents in interpreted R is orders of magnitude too slow. R's performance story is "push the
  loop into C/Fortran/vectorized primitives," not "compile your loop."

So a faithful port cannot be pure R. The performance-critical layer has to live in a
compiled language.

---

## 3. razer's answer: Rust + extendr + Rayon

`razer` fills both gaps with **Rust**, exposed to R through **extendr**:

| LASER (Python) | razer (R) | Role |
|---|---|---|
| NumPy array | **`Column`** — a Rust-owned `Vec<T>`, dtype-tagged (`i8/u8/i16/u16/i32/u32/f32/f64`), held in R as an opaque external pointer | typed, mutable, in-place agent/report storage |
| Numba `@njit(parallel=True)` kernels | **Rust kernels** parallelized with **Rayon**, compiled ahead-of-time at package install | the hot per-agent loops |
| `LaserFrame` (capacity + count + properties) | a `people` **environment** of `Column`s with `count`/`capacity` + `squash()` | the agent SoA container |
| node-report `LaserFrame` | 2-D (`n_ticks × n_nodes`) report `Column`s (`allocate_vector`) | per-node time series |

The defining trade-off of the `Column` approach: agent data is **not** an R-visible array.
It lives in Rust and crosses into R only on an explicit `$values()` snapshot. That buys
real wins — no copy-on-modify, narrow dtypes R can't express, zero-copy kernel access — at
the cost of the easy, interactive NumPy-style poking that LASER users enjoy.

Rust-vs-Numba is, on balance, a **favorable** swap: kernels are compiled once at install
(no per-session JIT warm-up), the borrow checker makes the in-place mutation safe, and
Rayon's work-stealing parallelism is mature. The cost is that **writing a new kernel
requires writing Rust** — there is no in-language escape hatch comparable to "just write
another `@njit` function in Python."

---

## 4. Component-by-component comparison

### 4.1 Data structures — *equivalent in spirit, different in ergonomics*

`LaserFrame` and the `Column`/`people`-environment pairing are the same idea: SoA, reserved
capacity, a live count, dead-slot reclamation via `squash`. Differences:

- **Visibility.** LASER's properties are NumPy arrays you can index and mutate anywhere.
  razer's are opaque handles; you read with `$values()`/`$col()` and write with
  `$set()`/`$set_col()`. razer is safer and more compact but less hackable at the REPL.
- **Dynamic property creation.** Both can add properties; LASER's `add_scalar_property` /
  `add_vector_property` on a frame is a touch more first-class than razer's "allocate a
  `Column` and stash it in the `people` env" idiom (which `run_model`'s `init` callback
  exists to host).
- **Persistence.** `laser-core` `LaserFrame` supports `save_snapshot` / `load_snapshot`.
  **razer has no snapshot mechanism** — a real gap for long runs and checkpoint/restart.

### 4.2 Parameters — *different*

LASER's `PropertySet` is a structured, serializable parameter bag passed to the model and
components. razer has **no `PropertySet` equivalent**; parameters are ordinary `run_model()`
arguments and values closed over in callbacks. Fine for the menagerie, weaker for
reproducible parameter sweeps and provenance.

### 4.3 Architecture — *the biggest divergence*

This is where the two part ways most.

- **laser-generic is a component framework.** A model is an ordered list of component
  objects; you extend it by writing a new component class and registering it. The set of
  expressible models is open-ended and the composition is explicit and introspectable.
- **razer is a fixed menagerie plus callbacks.** `run_model()` wires the eight closed-
  population models (SI/SEI/SIS/SEIS/SIR/SEIR/SIRS/SEIRS) in one correct per-tick order, and
  user extension happens through three lifecycle **callbacks** (`step_enter` / `step_update`
  / `step_exit`), `extra_states` (register a `"V"` or `"Q"` state), and the generic
  `step_timer_expire(_set)` kernels. This covers a lot — vaccination, waning, quarantine,
  importation, vital dynamics — and the worked examples prove it. But it is **composition by
  escape hatch, not by first-class component**. There is no plugin registry, no per-component
  init/report contract, no way for a third party to ship "a razer component" the way they can
  ship a laser-generic one. And because new *fast* dynamics need Rust, the callback route is
  capped at whatever you can do in R between kernel calls.

### 4.4 The per-tick model loop and the `R0 = β·D` insight — *razer is arguably more correct*

razer settled on a single per-tick ordering for **every** model —
`carry_forward → step → calc_foi → transmission` — with `calc_foi` reading the **settled
start-of-interval** infectious census `I[t]` (not the working column being built). Because
the FOI no longer depends on where the step kernel runs, every model uses one ordering and
an infectious agent contributes to the force of infection on exactly the `D` census columns
it occupies, giving the **full `R0 = β·D`** (never `β·(D−1)`) for both direct-S→I and
E-entry families, validated against the Kermack–McKendrick final size. This repository's
history shows the off-by-one (`β·(D−1)`) being found and fixed precisely because the SoA
ordering made the tally/infection interaction explicit. **Whether laser-generic realizes
`β·D` or `β·(D−1)` under its component ordering is worth auditing** — see §6.

### 4.5 Distributions — *equivalent, smaller surface*

razer ships `dist_normal/constant/uniform/gamma/poisson/beta/exp/lognormal/logistic` as
Rust-backed `Distribution` handles with thread-local sampling, validated against R's
`q*`/`d*` references. This mirrors LASER's use of parameterized distributions for timers.
The set is smaller and the `dist_` prefix exists only to dodge R's namespace collisions
(`stats::poisson`); functionally comparable.

### 4.6 Spatial coupling / migration — *at parity*

razer ports the `laser-core` migration family essentially 1:1: `distances` (haversine),
`gravity`, `radiation`, `stouffer`, `competing_destinations`, and `row_normalizer`, with the
network feeding `calc_foi`'s redistribution of the force of infection. The repository notes
these match laser-core's reference outputs to floating-point precision. This is the
strongest area of parity.

### 4.7 Demographics & initialization — *at parity, with one addition*

razer ports the alias-method age sampler (`AliasedDistribution` / `sample_pyramid_ages` /
`load_pyramid_csv`) and the `KaplanMeierEstimator` (date-of-death conditioned on current
age), plus `calc_capacity`. It **adds `calc_capacity_cdr`**, a mortality-aware bound on the
*peak living* population (vs. the cumulative-births bound), so a `squash`-reclaimed century
run can be sized realistically. That addition is a backport candidate (§6).

### 4.8 RNG & reproducibility — *razer is better*

razer's RNG (`set_seed` / `unset_seed`, `SmallRng`/xoshiro256++) is **reproducible
independent of CPU/thread count**: parallel kernels split agents into fixed-size chunks and
seed each chunk deterministically from `(call_base, chunk_index)`, so results don't depend on
how Rayon schedules work. Reproducibility that survives a change in core count is hard to get
from a Numba `prange` + global NumPy RNG design, and is a genuine improvement (§6).

### 4.9 Performance — *comparable, different cost curve*

Both end up running native, parallel kernels over contiguous typed arrays. Expected
differences: razer pays **zero JIT warm-up** (Rust is compiled at install; Numba compiles on
first call each session) but pays a **one-time Cargo build** at install. Steady-state
throughput should be in the same ballpark; the repository reports tens of millions of agents
processed per kernel in tens of milliseconds. No head-to-head benchmark exists in this tree,
so treat "comparable" as a design-level expectation, not a measured result.

---

## 5. What works well, and what doesn't

### Works well
- **The SoA port is faithful and the Rust substitution is sound.** `Column` + Rayon is a
  legitimate stand-in for NumPy + Numba, and arguably cleaner (compiled once, memory-safe,
  narrow dtypes).
- **Spatial, demographic, and distribution layers are at or near parity** and validated.
- **`run_model()` is a genuinely nice on-ramp.** One call for any menagerie model, with the
  correct ordering baked in, is more foolproof than hand-assembling components.
- **Reproducibility and the `β·D` ordering are improvements over the source design.**
- **Documentation is strong** — eight teaching articles, a reference site, validated examples.

### Doesn't work as well
- **No component framework.** Extensibility is callbacks + Rust, not first-class plugins. The
  open-ended composability that is laser-generic's whole point is absent.
- **Opaque data.** Agent arrays aren't R-visible; the interactive, "just `np.where` it"
  ergonomics are lost. R users expecting data-frame-like fluency will feel the friction.
- **Extending fast paths requires Rust.** A pure-R modeler hits a wall the moment they need a
  novel per-agent kernel; there is no in-language JIT to fall back on.
- **Missing infrastructure:** no `PropertySet`, no snapshots, no `SortedQueue`/event queue, no
  within-host or age-structured-contact machinery, a fixed model set.
- **Maturity and bus factor.** razer is `0.0.0.9000`, single-author, pre-release. LASER is an
  established, multi-contributor, actively developed toolkit.

---

## 6. Changes to give razer parity with LASER

Roughly in priority order:

1. **A component abstraction.** A registry of process objects (each with `init`/`step`/report
   hooks) over the `model` environment would convert the callback escape hatch into a real,
   third-party-extensible framework — the single biggest gap.
2. **A `PropertySet` equivalent** — a structured, serializable parameter object for provenance
   and sweeps.
3. **Snapshot save/load** for checkpoint/restart on long runs.
4. **An R-friendly view layer** over `Column`s (lazy accessors, `as.data.frame`, maybe an
   ALTREP-backed read-only view) to recover NumPy-style ergonomics without copies.
5. **More dynamics**: within-host, age-structured contact matrices, broader transmission
   options, an event/queue mechanism (`SortedQueue`).
6. **A documented Rust extension path** (or a constrained DSL) so adding a kernel is a
   supported workflow, not a fork.
7. **A real head-to-head benchmark** against LASER to replace design-level "comparable" claims
   with numbers.

## 7. Decisions in razer that should be backported to LASER

Stated plainly, these are places where the port improved on the original:

1. **Thread-count-independent reproducibility.** Seeding fixed-size agent chunks from
   `(call_base, chunk_index)` instead of relying on a global RNG under `prange` gives runs that
   reproduce regardless of core count. This is a correctness/credibility win any ABM should
   want, and LASER should adopt the pattern.
2. **`calc_capacity_cdr` (peak-living capacity bound).** If `laser-core` sizes only by
   cumulative births, the mortality-aware peak-living bound is a strictly better allocation for
   `squash`-reclaimed long runs.
3. **The `R0 = β·D` ordering discipline.** The "read the settled `I[t]`, place `calc_foi`
   immediately before transmission" rule eliminates a subtle `β·(D−1)` bias. laser-generic's
   component ordering should be audited against the Kermack–McKendrick final size; if it
   exhibits the same off-by-one, adopt the ordering fix.
4. **A single, documented per-tick ordering** that works across the whole family (rather than
   per-model sequencing) reduces a class of bugs and is worth importing as a convention.

---

## 8. Honest assessment: `razer` vs. LASER-through-`reticulate`

The blunt question for an R user: *should I use `razer`, or call LASER from R via
`reticulate`?*

**The case for `reticulate` + LASER (today, for serious work): stronger.** You get the
*actual* toolkit — the full component ecosystem, the breadth of models, the active
development and bug-fixing, NumPy arrays that `reticulate` will hand back to R as native
vectors/matrices, and parameters/objects you can introspect. The price is real but bounded:
a managed Python environment, marshaling across the R↔Python boundary, two-language
debugging, and `reticulate`'s occasional friction. For anyone doing production spatial ABM
of the kind LASER already supports, that price buys vastly more capability than `razer`
currently offers. **If your requirement is "do the modeling LASER does," `reticulate` + LASER
is the more capable answer right now, and it is not close on feature breadth.**

**The case for `razer`: narrower but real.** It wins decisively when:
- **A Python runtime is unacceptable** — locked-down environments, pure-R deployment targets,
  CRAN-style distribution, or teams that simply will not maintain a Python env.
- **The model is in the menagerie** (SI…SEIRS, ± vital dynamics, importation, simple
  interventions) — here `razer` is *easier* than wiring LASER components, fully native, and
  fast.
- **Teaching and onboarding** — the articles + `run_model()` make razer an excellent vehicle
  for explaining LASER's ideas to R-literate epidemiologists without a Python detour.
- **You value the specific improvements** razer made (reproducibility, the `β·D` discipline).

**The honest verdict.** `razer` is a *successful proof that the LASER architecture ports to
R*, and in a few spots it is cleaner than the original. It is **not** a substitute for LASER's
breadth, and pretending otherwise would do users a disservice. For most R users doing
non-trivial LASER-style work today, `reticulate` + LASER is the pragmatic choice — you get the
real, maintained toolkit. `razer` earns its place as (1) a Python-free option for the model
menagerie it covers, (2) a first-class teaching tool, and (3) a design laboratory whose better
ideas should flow back upstream. Its long-term relevance depends almost entirely on two
things: whether it grows a genuine **component framework** (without which it stays a fixed
menagerie), and whether enough R users actually need native-R ABM to justify maintaining a
parallel implementation rather than investing that effort in making `reticulate` + LASER
seamless. If neither holds, razer's most durable contribution may be the upstream
backports in §7 — and that would still make the experiment worthwhile.
