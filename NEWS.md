# razer (development version)

All notable changes to this project are documented here. (Headed as the development
version so pkgdown renders this file as the website's Changelog page; sub-sections follow
the "Keep a Changelog" convention.)

## Unreleased

### Changed
- **Terminology: "state", not "compartment".** razer is agent-based — each individual
  carries a *state* (susceptible, exposed, infectious, recovered, vaccinated, quarantined,
  maternally immune, deceased) — not compartmental (cohorts tracked only by a count). All
  prose, comments, docstrings, plot labels, and user-facing messages now say "state"
  accordingly. Two `run_model()` warning strings changed text: "model X has no E/R
  **compartment**" → "no E/R **state**", and the `extra_states` error now reads "must not
  repeat the model's own **states**" (matching test patterns updated).

### Examples / docs
- **Examples now play well in RStudio.** Every example script is device-aware — under
  `Rscript` it writes its PNG to `examples/output/` as before, but `source()`d in an
  interactive RStudio session it draws to the Plots pane instead (a small
  `to_png`/`open_png`/`close_png` helper). Added `examples/notebooks/getting_started.Rmd`,
  an R Notebook (R's Jupyter analog) walking through `run_model()`, a model comparison, and
  a callback-based intervention with inline plots; and an "Running in RStudio" section in
  the examples README. Added three more annotated teaching notebooks pairing with the
  scripts: `attack_fraction.Rmd` (R0 = β·D and the Kermack–McKendrick final size),
  `endemic_dynamics.Rmd` (the S*/N = 1/R0 equilibrium, importation, seasonality), and
  `interventions.Rmd` (custom states + callbacks: vaccination ± waning, quarantine).
  Completed the set with `spatial_metapopulation.Rmd` (migration networks + an England &
  Wales attack-rate map), `demographics.Rmd` (age pyramids + Kaplan–Meier dates of death),
  `vital_dynamics_measles.Rmd` (births/mortality + a maternal `M` state), and
  `long_runs_and_memory.Rmd` (`calc_capacity` vs `calc_capacity_cdr` + `squash()`), so every
  example now has a paired teaching notebook.
- **Richer example visuals:** `simple_sir.R` (previously plotless) now plots the national
  S/I/R trajectory and a geographic England & Wales attack-rate map; `quarantine.R` gained a
  2×2 figure (infectious baseline-vs-quarantine, cumulative cases averted, the Q
  state, and per-test detections).
- **`quarantine.R` enhancements:** the testing programme now runs on a configurable schedule
  (`test_period`, default fortnightly) and is **leaky** — each infectious agent is detected
  with probability `sensitivity` (default 0.80; 20% false negatives, no false positives).
- **Teaching notebooks are now pkgdown articles.** The eight notebooks moved from
  `examples/notebooks/*.Rmd` to `vignettes/articles/*.Rmd` and are rendered into the package
  website (the *Articles* menu, grouped *Foundations* / *Building richer models* /
  *Interventions & scale*) at <https://clorton.github.io/razer/>. They are website-only
  (`.Rbuildignore`d, so they don't run during `R CMD check`); the previously committed
  static `.html` renders were removed in favour of the site. `_pkgdown.yml` also gained the
  missing `step_timer_expire` / `step_timer_expire_set` reference entries.
- **`CHANGELOG.md` renamed to `NEWS.md`** (the conventional R name) so pkgdown renders it as
  the site's Changelog page; headed as the development version for pkgdown's news parser.

### Infrastructure
- **GitHub Actions pinned by commit SHA.** Both workflows (`pkgdown.yaml` and
  `R-CMD-check.yaml`) now reference every action by its full commit hash (with the version
  tag in a trailing comment) instead of a mutable tag, for supply-chain integrity;
  `pkgdown.yaml`'s `actions/checkout` bumped to v6 to match `R-CMD-check.yaml`.

### Added
- **`run_model()` gains `capacity` and `extra_states`** so the callbacks can express
  vital-dynamics and extra-state models, not just the closed-population menagerie:
  - `capacity` (default = initial population) reserves agent-array slots a `step_*`
    callback can activate with `births` / `import_infections` — i.e. lets the population
    grow (size it with `calc_capacity()` / `calc_capacity_cdr()`).
  - `extra_states` registers agent states beyond the model's own S/E/I/R: each gets a
    census Column that is carried forward, totalled into `N`, and seeded at tick 0 from
    agent states. A known [laser_states()] name (`"M"`) keeps its code; a NEW name (e.g.
    `"V"` for vaccinated) is assigned a free code, becoming a genuine new agent state, with
    the codes exposed in `model$states`. The disease kernels branch only on S/E/I/R/M, so
    an agent in a new state is left untouched (not infected, not transitioned) — its
    transitions are user-driven via callbacks (move `S`→`V` to vaccinate; for waning, set a
    timer and run `step_timer_expire(V, S)` in `step_update`), needing no kernel change.
    For `"M"`, run_model additionally applies the kernels' built-in M→S waning each tick
    (the `waning_m` flow), closing the desync that discarding `waned` would cause.
  These unblock expressing constant-population vitals, importation, growth, and a maternal
  `M` state through `run_model()` + callbacks (see the converted `simple_sir.R`,
  `endemic_sir*.R`, `long_run_squash.R`).

### Fixed (red-team follow-up)
- **`run_model` input validation & honesty.** Rejects a non-finite/negative `r0`; rejects
  `NA`/`Inf`/fractional/negative `population` and `E`/`I`/`R` seed columns with clear
  messages (no more silent truncation or cryptic base-R errors); and **warns** when a
  scenario seed column or a period argument is supplied that the chosen model does not use
  (e.g. an `E` column or `incubation_period` passed to `SIR`).
- **`step_update` docstring corrected:** it runs after the step kernel but before
  `calc_foi`, which reads the settled slice `t` — so census edits there affect the NEXT
  tick's FOI; edit the drivers (`beta`/`seasonality`/`network`) to affect the current tick.
- **`bincount_where()` / `bincount_where_wt()`: `count` is now required** (no longer
  defaults to the Column's capacity), so an over-allocated population can't silently tally
  reserved, inactive slots. Pass the active count (e.g. `people$count`).

### Changed
- **One per-tick ordering for every model — `calc_foi` immediately before `transmission`.**
  `calc_foi` now reads the **settled start-of-interval** infectious census `I[t]` instead
  of the working column `I[t+1]` (`src/rust/src/transmission.rs`). Because its value no
  longer depends on where the step kernel runs, every model uses the single ordering
  `carry_forward → step → calc_foi → transmission`, with no step kernel between `calc_foi`
  and `transmission`. An infectious agent still contributes to the FOI on exactly the `D`
  census columns it occupies, so `R0 = beta * D` holds for both direct-S→I and E-entry
  families with no per-family special-casing (the old dual ordering is gone). Re-validated
  against Kermack–McKendrick (SIR and SEIR attack fractions still match). `run_model` and
  all examples/tests updated to the single ordering.
- **`run_model()` now returns a `model` environment and takes lifecycle callbacks.** It
  bundles `$people`, `$nodes`, `$network`, and `$carry` into one environment (returned, and
  passed to the callbacks), and accepts optional `init(model)` (once, before the loop) plus
  three per-tick callbacks: `step_enter(model)` (start of tick), `step_update(model)`
  (between the step kernel and `calc_foi`), and `step_exit(model)` (end of tick). It now
  seeds **any**
  of the model's `E`/`I`/`R` states that the scenario supplies a column for (a state
  absent from the scenario, or from the model, is not seeded), and records the per-node
  flows `incidence` (all models), `onset` (E→I), `recovery` (I-exit), and `waning` (R→S).
  `$nodes$recoveries` is renamed `$nodes$recovery`.
- **Binning functions renamed for consistency:** `bincountw` → `bincount_wt`,
  `count_by_where` → `bincount_where`.

### Added
- **`bincount_where_wt()` — weighted, predicate-filtered binning by group**
  (`src/rust/src/bincount.rs`, `R/bincount.R`). The weighted twin of `bincount_where()`:
  sums a per-agent weight per group over the first `count` agents matching `prop <op>
  value`, in one parallel pass. Completes the family `bincount` / `bincount_wt` /
  `bincount_where` / `bincount_where_wt`. Covered by `tests/testthat/test-bincount-where.R`.
- **`calc_capacity_cdr()` — squash-aware capacity estimate** (`R/calc_capacity.R`). The
  mortality-aware companion to `calc_capacity()`: bounds the **peak living population** (net
  births minus a conservatively underestimated death rate) rather than the cumulative
  number ever born, so a [squash()]-reclaimed run can model decades/centuries without one
  slot per agent ever born. The death rate is underestimated by the safety factor (credited
  at `1/(1+safety_factor)`). Covered by `tests/testthat/test-calc-capacity-cdr.R` (1e6 over
  100 years, CBR 30 / CDR 15).
- **`bincount_where()` — flexible predicate-filtered agent binning by group**
  (`src/rust/src/bincount.rs`, `R/bincount.R`). A count-aware, filtered `bincount`: for
  each group `g` in `0..n_groups`, counts how many of the first `count` agents both have
  `group[i] == g` and satisfy `prop[i] <op> value` (`op` ∈ eq/ne/lt/le/gt/ge), in one
  parallel pass with no copy of the property into R. Answers queries like "exposed by
  node" (`prop = state, op = "eq", value = E`) or "under-fives by node"
  (`prop = dob, op = "gt", value = tick - 5*365`). Returns an integer vector for ad-hoc
  use, or writes into a report Column slice for per-tick model loops. Scans only the live
  prefix `count`, so it stays correct for variable-population (post-`squash`/birth) models.
  Covered by `tests/testthat/test-bincount-where.R`.
- **Agent compaction: `Column$squash(keep)` + the `squash(people)` helper**
  (`src/rust/src/column.rs`, `R/squash.R`). `Column$squash(keep)` stably compacts a 1-D
  Column to the elements flagged by a logical mask, returning the kept count;
  `squash(people)` applies one mask (default: still-alive, `state != D`) across every
  per-agent Column in a people environment so they stay row-aligned, and updates the
  active count. This reclaims the slots of agents that have left the simulation (e.g. the
  deceased) so the per-tick kernels stop iterating over them and `births` can refill the
  slots — the compaction the Column world had been missing. The per-node census is
  aggregate and unaffected. Covered by `tests/testthat/test-squash.R`.
- **Generic composable transition kernels** `step_timer_expire(from, to)` and
  `step_timer_expire_set(from, to, duration)` (`src/rust/src/steps.rs`). One timed
  `from_state → to_state` transition each (untimed vs timed destination), returning
  per-node counts — the building blocks for composing models beyond the named menagerie
  (a vaccinated `V`, a second infectious stage, …) from R without writing Rust. Covered
  by `tests/testthat/test-timer-expire.R`, including a hand-built SIRS assembled purely
  from the generics that conserves population. (Also confirmed `bincount` / `bincountw`
  remain Column-native and usable with the Column model — they tally a per-agent property
  Column into a per-node census/report slot; a model-use test was added.)
- **`run_model()` — a high-level runner for the closed-population menagerie**
  (`R/run_model.R`). One call builds the agent population and per-node census from a
  scenario, seeds it, and runs any of the eight SI/SEI/SIS/SEIS/SIR/SEIR/SIRS/SEIRS
  models on the Column kernels in the correct order — the `calc_foi` placement that
  yields `R0 = beta · D`, and every `move_count` census delta — so a model is one call
  instead of ~100 lines of hand-wiring, and the census cannot silently desync from the
  agents. Optional spatial `network`, `seasonality`, and a reproducible `seed`. Returns
  the `people`/`nodes` environments (census + `incidence`/`recoveries` flows). Validated
  against Kermack–McKendrick (SIR/SEIR attack fractions match to ~1e-3); covered by
  `tests/testthat/test-run-model.R` (all eight models, population conservation,
  reproducibility, input validation). (Hand-wired examples remain for models beyond the
  closed-population menagerie — vital dynamics, importation, maternal immunity.)
- **Reproducible, seedable RNG (`set_seed()` / `unset_seed()`).** All kernel randomness
  now flows through `src/rust/src/rng.rs`: after `set_seed(s)` an entire razer run is a
  deterministic function of `s` and the order of kernel calls, **independent of CPU/thread
  count**. The parallel kernels split agents into FIXED-size chunks (not one-per-thread)
  and seed each chunk deterministically from `(call_base, chunk_index)`, so the result is
  reproducible regardless of how Rayon schedules the work; a per-call counter gives each
  kernel invocation (and tick) an independent stream. Without a seed the behaviour is
  unchanged (entropy-seeded, random). The model RNG is `SmallRng` (xoshiro256++). All
  agent-loop kernels and the `Distribution` / `AliasedDistribution` / `KaplanMeierEstimator`
  samplers were converted off `thread_rng`. Covered by `tests/testthat/test-rng.R`.

### Changed
- **Unified the per-tick kernels into a return-counts menagerie (all `u16` timers).**
  The transmission and step kernels no longer take node census/flow buffers: they mutate
  the per-agent arrays and **return per-node counts**, which the model applies to the
  states it maintains via the new `move_count(from, to, counts, tick)` helper (so a
  model allocates only the states it has). The agent `timer` is now `u16` everywhere.
  - `transmission(state, timer, nodeid, count, foi, tick, to_state, duration)` →
    per-node infection counts (S→E or S→I, sets a u16 timer). Replaces the old
    `transmission` / `transmission_u16` (their `s_count`/`to_count`/`incidence` args are
    gone). New `transmission_si(state, nodeid, count, foi, tick)` → counts is the SI
    model's S→I into an absorbing `I` (no timer).
  - Three step kernels cover all of SI/SEI/SIS/SEIS/SIR/SEIR/SIRS/SEIRS, replacing
    `sir_step` and `measles_step`: `step_si` (M→S, E→I), `step_sir` (+ I→`absorbing_state`,
    S or R), `step_sirs` (+ I→R with immunity timer, R→S). Each is a single u16-timer pass
    leading with M→S and returns a named list of per-node transition counts.
  - `mortality(...)` → `list(m, s, e, i, r)` of per-node deaths by source state, and
    `births(...)` → `list(count, born)` (new active count + per-node births), instead of
    writing the census/flow directly. (`calc_foi` is unchanged; `constant_pop_vitals_sir`
    and `import_infections` keep writing their own census, now u16-timer-aware.)
  - All examples and tests migrated; the attack-fraction pair still validates
    Kermack–McKendrick with `R0 = beta · D`. Docs (`CLAUDE.md`, READMEs, `_pkgdown.yml`)
    updated with the menagerie table and the return-counts convention.

### Removed
- **The legacy `LaserFrame` world.** The package now has a single agent-storage and
  kernel stack — the `Column` typed arrays plus the per-tick Column kernels. Removed:
  the `LaserFrame` struct (`src/rust/src/laser_frame.rs`, `R/laser_frame.R`), the
  monolithic `epidemic.rs` step kernels (`step_transmission_si`/`_se`,
  `step_timer_expire`/`_set`, `step_exposed_ei`, `step_infectious_ir`/`_is`,
  `step_recovered_rs`, `step_mortality_cdr`, `step_births_cbr`), the `run_sir()` runner
  (`R/models.R`) and the `model-composition` help topic. `epidemic.rs` now provides only
  the shared state codes (`laser_states()`). The monolithic `step_transmission_si` was
  the one kernel that could not realize the full `beta · D` (see Fixed); retiring it
  removes the last `beta · (D − 1)` path. `examples/sir_attack_fraction.R` and
  `seir_attack_fraction.R` were rewritten on the Column kernels (both now validate
  Kermack–McKendrick with `R0 = beta · D`); SEIR reuses `measles_step` with the maternal
  state empty. Associated tests removed; `laser_states` keeps its own test.

### Fixed
- **SIR models now realize the full `R0 = beta · D`, not `beta · (D − 1)`.** With direct
  `S→I` transmission a newly infectious agent enters `I` *after* the force-of-infection
  tally; if recoveries were also excluded from the tally the agent contributed on only
  `D − 1` ticks. In the Column-kernel examples (`simple_sir.R`, `endemic_sir.R`,
  `endemic_sir_seasonal.R`) the per-tick order is changed to run **`calc_foi` before
  `sir_step`** (recovery), so an agent is counted on its recovery tick and the effective
  infectious period is the full `D`. The endemic SIR now settles with `S` at `1/R0`
  (was ~`1/(beta·(D−1))`). `calc_foi`'s docstring and the `CLAUDE.md` modeling notes are
  updated to document the ordering and to forbid `beta·(D−1)` models. (The legacy
  monolithic LaserFrame `step_transmission_si`, used by `run_sir()` and
  `sir_attack_fraction.R`, still yields `beta·(D−1)` — a structural limitation of
  computing the tally and infection in one kernel; those are slated to move to the
  Column kernels.)

### Changed
- **Breaking:** `step_transmission_si()` and `step_transmission_se()` now take a
  required fifth argument, `network`, enabling spatial coupling of contagion (see
  Added). razer is spatial-first: there is no "no coupling" special case — pass an
  **all-zero matrix** to run nodes independently (identical to the previous
  per-node force of infection). All existing callers were updated.

### Added
- `transmission_u16(state, timer, nodeid, count, foi, S, to_count, incidence, tick, to_state, duration)`
  (`src/rust/src/measles.rs`) — the measles S→E transmission step: identical to
  `transmission()` but writes a **uint16** timer (the incubation clock). Reuses the
  per-node `1 - exp(-foi)` probability and the parallel-tally pattern. Covered by
  `tests/testthat/test-transmission-u16.R` (incl. a 1M-agent parallel-vs-census check).
- `births(state, timer, nodeid, dob, dod, count, birth_rate, M, births_flow, maternal_duration, km, tick)`
  (`src/rust/src/births.rs`) — crude-birth-rate births. Each living agent births with
  per-node probability `1 - exp(-birth_rate)`; each newborn activates a reserved slot as
  `M` (maternal immunity) with a `timer` from `maternal_duration`, `dob = tick`, and a
  `dod` drawn from the `KaplanMeierEstimator` (`tick` + a newborn age at death). Returns
  the grown active count; updates the `M` census at `tick+1` and the `births` flow at
  `tick`; caps at capacity. The `KaplanMeierEstimator` gained a crate-internal
  `sample_newborn_age_at_death`, and `Column` gained `as_u8` / `as_u32_mut`. Covered by
  `tests/testthat/test-births.R`.
- `examples/engwal_measles.R` is now a **complete single-patch measles model**: the
  per-tick loop wires `carry_forward_states → measles_step → calc_foi →
  transmission_u16 → births → mortality`, with initial infectious timers seeded and a
  CBR-driven `birth_rate` grid. Over 20 years it reproduces the classic measles
  inter-epidemic cycles (damping toward the endemic equilibrium) with susceptibles
  regulated at the `N/R0` threshold; writes a reservoirs/incidence dynamics plot.
- `measles_step(state, timer, nodeid, count, M, S, E, I, R, inf_duration, tick)`
  (`src/rust/src/measles.rs`) — the combined timed-transition kernel for the measles
  model, advancing all three timed states in a single pass over a **uint16**
  timer: **M→S** (maternal-immunity waning → Susceptible), **E→I** (incubation end →
  Infectious, drawing a fresh infectious-period timer from `inf_duration`), and **I→R**
  (recovery → Recovered, lifelong). Branching on each agent's entry state means every
  agent is touched at most once per tick, so a just-arrived `I` is not also recovered
  the same tick — the downstream-first ordering achieved structurally rather than by
  sequencing separate kernels. The decrement is guarded so the u16 timer never
  underflows. Census deltas are applied at column `tick+1`; parallelized across cores
  with private per-thread node buffers. Covered by `tests/testthat/test-measles-step.R`
  (incl. a parallel-vs-census check at 900k agents).
- **Natural (non-disease) mortality by date of death** for the measles model.
  - New `M` (maternal-immunity) state: `STATE_M = 4` and `laser_states()` now
    returns `c(S=0, E=1, I=2, R=3, M=4, D=-1)` (`src/rust/src/epidemic.rs`).
  - `mortality(state, dod, nodeid, count, M, S, E, I, R, deaths, tick)`
    (`src/rust/src/mortality.rs`) — retires every living agent whose date of death
    `dod` (an absolute `u32` tick) has arrived (`dod <= tick`): sets its state to `D`,
    decrements the M/S/E/I/R census it occupied at column `tick+1`, and adds to the
    per-node `deaths` flow at column `tick`. Parallelized across cores with private
    per-thread node buffers. Covered by `tests/testthat/test-mortality.R`.
  - `Column` gains a crate-internal `as_u32` accessor (`src/rust/src/column.rs`).
  - `examples/engwal_measles.R` now allocates initial ages from an age-distribution
    curve, draws each agent's age at death with a `KaplanMeierEstimator` built on the
    SAME curve, stores `dob` (signed `i32`, = −age) and `dod` (`u32`), and runs the
    `mortality()` kernel each tick. Writes an age-curve/life-table plot and a
    living-population/deaths plot to `examples/output/`.
- `calc_capacity(birthrates, initial_pop, safety_factor = 1)` (`R/calc_capacity.R`) —
  estimate the per-node agent capacity to preallocate for a population growing under a
  (possibly time-varying) crude birth rate, ported from laser.core. Treats births as
  geometric growth (daily rate `(1 + CBR/1000)^(1/365) - 1` summed over steps and
  exponentiated) with a `safety_factor` headroom term. Returns whole-valued R doubles
  (which represent integers exactly to `2^53`); unlike laser.core's `uint32` array it
  does **not** clamp, but `warning`s if any per-node estimate exceeds
  `.Machine$integer.max` — the `i32` `count` limit of razer's allocators. Accepts
  `birthrates` as an `nsteps x nnodes` matrix or a 2-D razer `Column` (e.g. from
  `values_map`). Covered by `tests/testthat/test-calc-capacity.R`.
- **Demographics: realistic age structure and dates of death**, ported from
  laser.core.
  - `AliasedDistribution` + `aliased_distribution(counts)` (`src/rust/src/pyramid.rs`)
    — a Vose alias-method sampler over 0-based bin indices weighted by integer counts
    (e.g. per-age-band populations). O(1) per draw, all-integer construction (no
    floating-point round-off). Methods `$sample_one()`, `$sample_n(n)` (parallelized
    across cores, each with a thread-local RNG), `$n_bins()`, `$total()`, `$alias()`,
    `$probs()`. Covered by `tests/testthat/test-aliased-distribution.R`.
  - `KaplanMeierEstimator` + `kaplan_meier_estimator(cumulative_deaths)`
    (`src/rust/src/kmestimator.rs`) — given a non-decreasing cumulative-deaths-by-year
    life table, samples a year or age (in days) of death conditioned on survival to a
    given age (Kaplan–Meier). Methods `$predict_year_of_death(ages_years, max_year)`,
    `$predict_age_at_death(ages_days, max_year)` (both parallelized), and
    `$cumulative_deaths()`. A negative `max_year` selects the last year in the table.
    Covered by `tests/testthat/test-kmestimator.R`.
  - R helpers (`R/pyramid.R`): `load_pyramid_csv(file)` parses a laser.core-format
    pyramid CSV (`"Age,M,F"` header, `"low-high,M,F"` bands, open-ended `"max+"`) into a
    `start`/`end`/`M`/`F` integer matrix; `sample_pyramid_ages(pyramid, n)` draws
    per-agent ages in days (band by population via the alias sampler, then a uniform day
    within the band). Covered by `tests/testthat/test-pyramid.R`.
  - `examples/aged_population.R` (+ `examples/data/pyramid_example.csv`) — end-to-end
    demo: load a pyramid, generate one million age-structured agents, build a synthetic
    Gompertz life table, and assign each agent a Kaplan–Meier date of death conditioned
    on its current age. Writes a two-panel age-structure / age-at-death plot to
    `examples/output/aged_population.png`.
- `import_infections(state, timer, nodeid, count, I, importations, sched_tick, sched_node, sched_count, duration, tick)`
  — schedule-driven importation of new infectious cases (`src/rust/src/vitals.rs`).
  For every schedule entry whose `sched_tick == tick` it activates `sched_count`
  RESERVED agent slots (those past the live `count`, in capacity-sized property
  arrays) in node `sched_node`: each is set Infectious with a `timer` drawn from the
  `duration` Distribution and its `nodeid` set. Per-node imports are added to the I
  census at column `tick+1` and written to the `importations` flow at `tick`; returns
  the grown live count (the caller stores it back into `people$count`). Asserts the
  imports fit the allocated capacity. Sequential (touches only the imported slots).
  Covered by `tests/testthat/test-import-infections.R`.
- `carry_forward_states(carry, tick, total = NULL, summands = carry)` — R convenience
  wrapper (`R/carry_states.R`). Carries each 2-D census Column in `carry` forward one
  tick (column `tick` -> `tick+1`, via `carry_forward`) and, if `total` is supplied,
  sets `total`'s column `tick+1` to the elementwise sum of the `summands` Columns at
  `tick+1` — e.g. carry S/I/R forward and total them into N (the population / FOI
  denominator) in one call, keeping N current as births, deaths, and imports change
  the states. Covered by `tests/testthat/test-carry-forward-states.R`.
- `$col(slot)` / `$set_col(slot, values)` accessor methods on `Column`
  (`src/rust/src/column.rs`) — read or overwrite one column (e.g. all nodes for one
  tick) of a 2-D Column in place, widening / casting like `$values()` / `$set()`.
  `carry_forward_states` builds on these. Covered by `tests/testthat/test-allocation.R`.
- `examples/endemic_sir.R` — an endemic two-patch SIR model (500,000 agents per patch,
  R0 = 3, seeded at the endemic susceptible fraction 1/R0 with the rest immune, a high
  crude death rate of 20 with no spatial/temporal variation, a small 1% inter-patch
  coupling) run for 10 years. Periodic importations spark transmission so a stochastic
  fade-out doesn't end the run; vital turnover resupplies susceptibles. Demonstrates
  `import_infections`, `carry_forward_states` (S/I/R carried and totalled into N each
  tick), and `constant_pop_vitals_sir` wired together; the susceptible fraction settles
  near 1/R0. The `run_endemic_sir()` runner takes a `progress` flag (text progress bar
  over the per-tick loop) and returns the `people`/`nodes` environments; the script
  then writes an S/I/R-over-time plot to `examples/output/endemic_sir_SIR.png` (the
  shared example output directory).
- `examples/endemic_sir_seasonal.R` — the endemic two-patch SIR model with a gentle
  annual sinusoid on the transmission coefficient (`beta * (1 + A*sin(2*pi*day/365))`).
  An amplitude sweep (A = 0.05–0.25, two stochastic runs each at 500k agents/patch)
  showed that seeding at the endemic susceptible fraction S ≈ N/R0 leaves R_eff ≈ 1
  (critically poised), so even a small forcing is amplified into large, phase-locked
  annual epidemic waves (peaks recur near the same day each year; autocorrelation at
  lag 365 ≈ 0.5). **A = 0.10 (±10%)** is the chosen default: the smallest forcing that
  gives clear annual waves while the per-year prevalence trough never collapses to zero
  — A ≥ 0.15 occasionally faded a patch out in the off-season (the importation floor
  then re-sparks it). The `run_endemic_sir()` runner gains a `seasonality` argument (any
  `values_map`-broadcastable shape); the script writes a two-panel forcing-vs-prevalence
  plot to `examples/output/endemic_sir_seasonal.png`.
- `constant_pop_vitals_sir(state, timer, nodeid, count, rate, S, I, R, births, deaths, tick)`
  — constant-population SIR vital dynamics (`src/rust/src/vitals.rs`). Each agent
  dies with probability `1 - exp(-rate[node])` (the caller passes a per-node daily
  death HAZARD rate grid — e.g. crude death rate `cdr / 1000 / 365` via `values_map`)
  and is immediately reborn susceptible (`state -> S`, `timer -> 0`); every event is
  recorded as both a birth and a death (equal under constant population) in the
  `tick` column of the `births`/`deaths` flow reports. The S/I/R node census is kept
  exactly in sync at column `tick+1`: a death out of I/R moves that count to S (a
  death out of S nets to zero), so `S+I+R` stays equal to the population and matches
  a direct agent census. SIR-specific (knows the S/I/R states). Parallelized
  across cores with private per-thread node buffers reduced at the end. Covered by
  `tests/testthat/test-vitals.R` (incl. a parallel-vs-census check at 1M agents).
  `examples/simple_sir.R`'s `run_sir_model()` gains a `cdr` argument and runs it as
  the last per-tick step, resupplying susceptibles for endemic dynamics.
- **Spatial SIR transmission + an incrementally-maintained node census** over
  `Column` buffers (`src/rust/src/sir.rs`). All three kernels parallelize the
  per-agent work across cores with private per-thread node accumulators (no shared
  writes) reduced at the end:
  - `calc_foi(infected, population, beta, seasonality, network, foi, tick)` —
    per-node force of infection. Frequency-dependent local rate
    `r[k] = beta[k]·seasonality[k]·I[k]/N[k]`, redistributed through the network as
    `r[k]·(1−ΣW[k,·]) + Σr[i]·W[i,k]`, written to `foi[tick]`. Reads the
    **post-recovery** census `I[tick+1]` and `N[tick+1]` (the working columns after
    `carry_forward` + `sir_step`, so this tick's recoveries are excluded and the
    denominator is the current — eventually vital-dynamics-varying — population),
    and the exogenous `beta`/`seasonality` modifier grids at the interval column
    `tick`. `beta` and `seasonality` are `n_ticks × n_nodes` grids built by
    `values_map` (below). Covered by `tests/testthat/test-calc-foi.R`.
  - `sir_step(state, timer, nodeid, count, I, R, recoveries, tick)` — recovery.
    Recovers expired infectious agents and applies the I→R delta to column `tick+1`
    of the census (and the per-node flow to `recoveries[tick]`). The caller seeds
    column `tick+1` first via `carry_forward` (below), keeping the census
    incremental (`count[t+1] = count[t] ± delta`) with no per-step re-census.
  - `transmission(state, timer, nodeid, count, foi, S, to_count, incidence, tick, to_state, duration)`
    — infection, generalized over the receiving state. Converts `foi[tick]` to a
    per-node probability `1−exp(−foi)` **once per node** (≈3× faster than recomputing
    `exp` per susceptible — benchmarked at 30.18M agents / 954 nodes / ~95%
    susceptible), and moves susceptibles into `to_state` (`I` for SIR or `E` for
    SEIR), drawing each agent's `timer` for that state from `duration` (incubation
    for `E`, infectious period for `I`). Applies the S→`to_state` delta to column
    `tick+1` of the `S` and `to_count` census (flow to `incidence[tick]`). The caller
    passes the matching `to_count` census, `to_state` code, and `duration` — and the
    `timer` column for whichever clock the receiving state uses, so incubation /
    infectious / waning durations can be tracked in separate timer vectors.
  Each kernel has a test asserting the parallel per-node accumulation matches a
  serial census of the resulting agent states at 1M agents / 200 nodes
  (`test-sir-step.R`, `test-transmission.R`). `examples/simple_sir.R`'s
  `run_sir_model()` gains `inf_duration`, `r0` (beta = `r0 / mean(inf_duration)`),
  and `seasonality` arguments, builds per-node population (`nodes$N`) and `beta` /
  `seasonality` modifier grids via `values_map`, seeds the census + agents from
  optional `scenario$I` / `scenario$R` columns, and runs the
  downstream-first per-tick loop `carry_forward → sir_step → calc_foi →
  transmission`. `nticks` is the number of recorded states (0..nticks-1; state 0 is
  the seeded initial, state nticks-1 the final), so dynamics run **nticks-1 times**
  (the intervals between states — no step on the final state, which would run off
  the window); the S/I/R census is `nticks × n_nodes` and the foi/recoveries/
  incidence flows are `(nticks-1) × n_nodes`. Verified that `S+I+R == population`
  per node at every state and the census matches a direct agent census (30.18M
  agents, ~0.5 s).
- `values_map(value, n_ticks, n_nodes)` (`R/grid.R`) — expands a flexible
  per-time and/or per-node `value` into a full `n_ticks × n_nodes` f64 `Column`: a
  scalar (constant), a length-`n_nodes` vector (per-node, constant over time), a
  length-`n_ticks` vector (per-tick, constant over space), or an `n_ticks × n_nodes`
  matrix (both). Used to build `calc_foi`'s `beta` and `seasonality` grids from
  whatever shape the modeler supplies. Covered by `tests/testthat/test-grid.R`.
- `carry_forward(counter, tick)` — copies column `tick` of a 2-D report `Column`
  onto column `tick+1` (any dtype), seeding the next tick so a dynamics kernel can
  update it in place. Called once per census state that must persist across ticks —
  S, I, R for SIR; add E for SEIR; users can carry their own states (e.g. a "V"
  vaccinated count). Extracted from `sir_step` so the census carry-forward is
  model-agnostic. Covered by `tests/testthat/test-carry-forward.R`.
- `allocate_vector(dtype, n_slices, slice_len)` — allocates a 2-D `Column` report
  buffer of `n_slices × slice_len` elements, SLICE-MAJOR (row-major) so each slice
  is a contiguous run: slice `s` is `s*slice_len .. (s+1)*slice_len`. The
  conventional use is a time-series report with the **first dimension time, second
  node** (`allocate_vector(dtype, n_ticks, n_nodes)`), so each tick's per-node
  values are contiguous and a step kernel can fill one tick's slice with no
  gather. `$values()` reads it back as an `n_slices × slice_len` (e.g.
  `n_ticks × n_nodes`) R matrix (transposed on copy, since R is column-major), so
  row `t` is tick `t`'s vector. `Column` carries `n_slices`/`slice_len` shape;
  `allocate_scalar` is the `n_slices == 1` case. Covered by
  `tests/testthat/test-allocation.R`.
- `bincount(values, nbins, counts, slot = 0)` — a parallel, NumPy-`bincount`-style
  histogram over a [Column] (`src/rust/src/bincount.rs`). For each bin `b` in
  `0..nbins` it counts how many elements of the integer-typed `values` Column equal
  `b` and writes the result into **slice `slot`** of the caller-provided `counts`
  Column (the entries `slot*slice_len .. slot*slice_len + nbins`), zeroing then
  accumulating. For a scalar `counts` the only slice is `slot = 0` (the whole
  vector); for a 2-D report from [allocate_vector()] `slot` selects a tick's row —
  so a per-tick loop can tally straight into the report with no extra buffer. Each
  Rayon worker accumulates into a private per-thread histogram (no shared-bin write
  collisions). One generic kernel serves every value width/signedness (`i8`..`u32`)
  and any numeric `counts` type via traits. Counts 30.18M agents into 954 node bins
  in ~15 ms; covered by `tests/testthat/test-bincount.R` (incl. a
  parallel-vs-`tabulate` cross-check). `slot` defaults to 0 via a thin R wrapper
  (`R/bincount.R`) over the extendr kernel, since extendr can't emit R-side default
  arguments.
- `bincountw(values, weights, nbins, counts, slot = 0)` — the weighted counterpart:
  sums each element's weight into its bin (`counts[b] = Σ weights[i]` over `i` with
  `values[i] == b`), à la `numpy.bincount(values, weights=...)`. Same parallel,
  collision-free, zero-then-accumulate design and `slot` targeting; `weights` may be
  any numeric `Column` (signed, unsigned, or floating point — all widened to f64 for
  accumulation), and the value×weight type dispatch is macro-generated over one
  generic kernel so weights are read in place with no copy. Covered by
  `tests/testthat/test-bincount.R` (incl. a parallel-vs-serial weighted-tally
  cross-check).
- `Column` and `allocate_scalar(dtype, count)` — a Rust-owned, dtype-tagged 1-D
  property array exposed to R as an opaque external-pointer handle
  (`src/rust/src/column.rs`). `allocate_scalar()` returns a zero-filled `Column`
  of any of eight element types — `i8`, `u8`, `i16`, `u16`, `i32`, `u32`, `f32`,
  `f64` (with aliases like `"integer"` = i32, `"double"` = f64, `"uint8"` = u8) —
  none of which (beyond i32/f64/u8) R can represent natively. The data lives in a
  Rust `Vec<T>` so it is exactly `sizeof(T)` per element and the step kernels can
  borrow `&mut [T]` and mutate it in place with no R copy-on-modify. Methods:
  `$length()`, `$dtype()`, `$fill(value)`, `$set(values)`, and `$values()` (an
  on-demand snapshot copied into the nearest native R vector — `integer` for the
  narrow integer types, `double` for `u32`/`f32`/`f64`). Covered by
  `tests/testthat/test-allocation.R`. `examples/simple_sir.R`'s `run_sir_model()`
  builds a lightweight `people` environment (`count`, `capacity`, a `u8` `state`
  Column, and a `u16` `nodeid` Column initialized 0-based from the patch
  populations, sized to the total population). Node ids are 0-based (0..N-1) to
  match the Rust kernels' direct per-node indexing; R-side joins to `scenario`
  rows add 1.
- **Migration-network models** ported from laser-core (`src/rust/src/migration.rs`),
  each taking a population vector and a symmetric `N × N` distance matrix and
  returning an `N × N` migration-weight matrix (zero diagonal, generally
  asymmetric):
  - `gravity(pops, distances, k, a, b, c)` — gravity model
    `k · p_i^a · p_j^b / d_ij^c`.
  - `radiation(pops, distances, k, include_home)` — radiation model
    (Simini et al., Nature 2012).
  - `stouffer(pops, distances, k, a, b, include_home)` — Stouffer's intervening
    opportunities (1940).
  - `competing_destinations(pops, distances, k, a, b, c, delta)` — Fotheringham's
    competing-destinations adjustment to the gravity model (1984).
  - `row_normalizer(network, max_rowsum)` — proportionally caps each row sum at
    `max_rowsum`, e.g. to bound the exported fraction before using a matrix as the
    transmission `network`.
  Outputs match laser-core's reference implementation to floating-point precision;
  covered by `tests/testthat/test-migration.R`. `examples/simple_sir.R` uses
  `radiation()` + `row_normalizer()` to build its spatial coupling network.
- `distances()` — a port of laser-core's `distance` (all-pairs case): given
  vectors of latitudes and longitudes (decimal degrees), returns the symmetric
  `N × N` great-circle (haversine) distance matrix in kilometres, with a zero
  diagonal. Uses a 6371 km mean Earth radius to match laser-core, validates the
  coordinate ranges, and fills the matrix column-by-column across Rayon worker
  threads (`src/rust/src/migration.rs`). Covered by
  `tests/testthat/test-distances.R`.
- `examples/simple_sir.R` — a worked SIR example built on the high-level
  `run_sir()` runner. Its setup loads the England & Wales measles patches as the
  node scaffold and builds the pairwise `distances()` matrix from their
  latitude/longitude (the geographic input for the spatial coupling network);
  the SIR wiring follows.
- `examples/data/EnglandWalesMeasles_places.csv` and
  `examples/data/convert_measles.py` — a shareable, one-row-per-patch CSV (name,
  initial 1944 population, latitude, longitude) for the 954 England & Wales
  registration districts, plus the Python converter that flattens it from the
  source `EnglandWalesMeasles.py` dataset.
- **Spatial coupling of contagion** via the new `network` argument on
  `step_transmission_si()` and `step_transmission_se()`, porting laser-generic's
  force-of-infection redistribution model. The `network` is an `n_nodes × n_nodes`
  matrix whose `[i, j]` entry is the fraction of node *i*'s force of infection
  exported to node *j*; the per-node coupled rate is
  `lambda[k] = r[k]·(1 - rowSums(W)[k]) + (t(W) %*% r)[k]` with local rate
  `r[k] = beta·I[k]/N[k]`, converted to a probability `p[k] = 1 - exp(-lambda[k])`.
  Total force of infection is conserved and the diagonal cancels (self-export has
  no effect). An all-zero matrix leaves `lambda[k] = r[k]` (independent nodes) —
  a convenient "poor man's" batch of parallel single-node runs in one call.
  `run_sir()` gains a required `network` parameter that validates
  shape / non-negativity / off-diagonal row-sum ≤ 1 and threads the matrix
  through. Covered by `tests/testthat/test-transmission-network.R` (directional
  leak, unconnected-node isolation, FOI-magnitude and diagonal-cancellation checks
  against the formula, all-zero independence, shape validation) plus additional
  `run_sir()` network cases.
- `run_sir()` — a high-level **model runner** that assembles the per-tick step
  kernels and their parameters into a single call (`R/models.R`). It takes a node
  data.frame (one row per node, an integer `population` column), builds an agent
  `LaserFrame` sized to the total population, assigns each agent to its node,
  seeds initial infections, and runs the downstream-first SIR loop (recovery I→R
  before transmission S→I). It records per-node state trajectories (`S`,
  `I`, `R`) and per-tick flows (`incidence` S→I, `recovery` I→R) into a node-level
  report `LaserFrame` attached as `attr(model, "report")`, alongside `runtime`,
  `nticks`, `model`, and `parameters` attributes. `infectious_period` accepts a
  `Distribution` or a bare number (promoted to `dist_constant`), and `progress`
  draws a text progress bar. Intended as the reference template for `run_seir()`,
  `run_si()`, `run_sis()`, … Covered by `tests/testthat/test-run_sir.R`
  (node assignment, seeding, population conservation, flow-vs-delta identities,
  cross-node isolation, attribute surface, and argument validation).
- Two **generalized timer-expiry kernels** (mirroring laser-generic's
  `nb_timer_update` and `nb_timer_update_timer_set`), parameterized by the
  `from`/`to` state codes so any timer-driven transition can be expressed without
  a bespoke kernel:
  - `step_timer_expire(people, from_state, to_state)` — transition into an
    *absorbing* (untimed) state; the timer is left at 0 on arrival (e.g. I→S, R→S).
  - `step_timer_expire_set(people, from_state, to_state, duration_dist)` —
    transition into a state with *its own* duration; a fresh per-agent timer is
    drawn from `duration_dist` on arrival (e.g. E→I, I→R with waning).
  The four named kernels (`step_exposed_ei`, `step_infectious_ir`,
  `step_infectious_is`, `step_recovered_rs`) are now thin, fixed-state wrappers
  over these two helpers, so the shared decrement/transition loop lives in one
  place. Covered by `tests/testthat/test-timer-kernels.R`, which proves the
  wrappers are equivalent to the generalized kernels.
- `distributions` module exposing parameterized probability distributions to R as
  opaque `Distribution` handles. Constructors use a `dist_` prefix to avoid masking
  base/stats functions (e.g. `base::gamma`, `stats::poisson`):
  - `dist_normal(mean, variance)` — Gaussian (second argument is the variance σ²,
    not the standard deviation).
  - `dist_constant(value)` — degenerate distribution; a fixed-value drop-in.
  - `dist_uniform(low, high)` — continuous uniform on `[low, high)`.
  - `dist_gamma(shape, scale)` — gamma in the shape–scale (k, θ) parameterization
    (mean `shape*scale`, variance `shape*scale^2`); strictly positive draws.
  - `dist_poisson(lambda)` — Poisson with rate/mean `lambda`; non-negative integer
    counts.
  - `dist_beta(alpha, beta)` — beta on `(0, 1)` with shape parameters α, β.
  - `dist_exp(rate)` — exponential with rate λ (mean `1/rate`).
  - `dist_logistic(location, scale)` — logistic (mean `location`, variance
    `scale^2·π²/3`); sampled by inverse-CDF transform.
  - `dist_lognormal(meanlog, sdlog)` — log-normal with log-space parameters
    (matches R's `qlnorm`).
- `Distribution$sample_one()` and `Distribution$sample_n(n)` — draw one / a batch
  of floating-point samples with a thread-local RNG (for interactive use, testing,
  and statistical validation). All draws are doubles.
- Internal Rust sampler `Distribution::sample` (Pattern B): the caller supplies the
  RNG, so one `&Distribution` can be shared by reference and sampled concurrently
  across Rayon worker threads.
- `tests/testthat/test-distributions.R` covering the constructors and their use in
  the step kernels, and `tests/testthat/test-distributions-validation.R` validating
  parameter wiring against R's reference implementations (`qnorm`, `qunif`,
  `qgamma`, `qpois`/`dpois`) over one million draws.
- `examples/sir_attack_fraction.R` and `examples/seir_attack_fraction.R` —
  runnable SIR and SEIR examples that plot the state trajectories and
  compare the simulated final attack fraction against the Kermack–McKendrick
  final-size relation `A = 1 - exp(-R0 * A)` across an `R0` sweep, with timing
  output (with `examples/README.md` and sample output plots).
- `CLAUDE.md` documenting the downstream-first transition-ordering convention for
  composing models from the step kernels.

### Changed
- Model compositions (SIR/SEIR/SEIRS test helpers and the examples) now call the
  per-tick transitions **downstream-first** (out of each timed state before
  into it: R→S, I→R, E→I, S→E) so that exposed/infectious periods reflect their
  full configured duration rather than being shortened by one tick. See `CLAUDE.md`.
- Distributions return **floating-point** values; the integer rounding for state
  timers now lives in the step kernels (a shared `duration_ticks` helper rounds to
  the nearest tick and clamps to a minimum of 1). `Distribution` no longer exposes
  an integer `sample_duration`.
- Every step that transitions an agent **into** a timed state takes a `Distribution`
  for that state's duration instead of a fixed `i32`. The duration is drawn per
  agent at the moment of transition. Pass `dist_constant(d)` for fixed-duration
  behavior or e.g. `dist_normal(mean, variance)` for a stochastic period:
  - `step_transmission_si(people, nodes, beta, inf_dist)` — was `inf_duration`
    (sets the **I** timer on S→I).
  - `step_transmission_se(people, nodes, beta, exp_dist)` — was `exp_duration`
    (sets the **E** timer on S→E).
  - `step_exposed_ei(people, inf_dist)` — was `inf_duration`
    (sets the **I** timer on E→I).
  - `step_infectious_ir(people, imm_dist)` — was `imm_duration`
    (sets the **R** waning timer on I→R; use `dist_constant(0)` for SIR/SEIR where
    the R timer is never read).
  Steps transitioning into untimed states (`step_infectious_is`,
  `step_recovered_rs`, `step_mortality_cdr`, `step_births_cbr`) are unchanged.
  All SI/SIR/SEIR/SEIRS test helpers updated to wrap their durations in
  `dist_constant()`.
