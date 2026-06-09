# razer (development version)

Razer is in active pre-release development. This page summarizes the capabilities of the
current development version rather than every intermediate change; the detailed,
commit-by-commit history leading to the first release lives in the git log. (The file is
headed as the development version so pkgdown renders it as the website's Changelog page.)

## Agent and report storage

- **`Column`** — a Rust-owned, dtype-tagged array (`i8`/`u8`/`i16`/`u16`/`i32`/`u32`/`f32`/
  `f64`) that R holds as an opaque external-pointer handle, so the kernels borrow and mutate
  it in place with no copy. `allocate_scalar(dtype, count)` makes a 1-D per-agent array;
  `allocate_vector(dtype, n_ticks, n_nodes)` makes a 2-D, slice-major time × node report
  buffer (each tick's per-node row contiguous). Accessors `$values()` / `$set()` / `$col()` /
  `$set_col()` / `$length()` / `$dtype()` move data to and from the nearest native R vector.
  Each agent's disease `state` is a `u8` (`laser_states()` returns
  `c(S = 0, E = 1, I = 2, R = 3, M = 4, D = -1)`, with `D` stored as `255`) and its state
  `timer` is a `u16` everywhere (maternal / immunity periods exceed a `u8`'s 255).
- **Incremental census** — `carry_forward(column, tick)` copies a report Column's column
  `t` onto `t+1`; `carry_forward_states(carry, tick, total)` carries a set of census Columns
  forward and re-totals them into `N` in one call, so the census is maintained as
  `count[t+1] = count[t] ± deltas` with no per-tick re-census.

## Distributions

- **`Distribution`** handles built by `dist_`-prefixed constructors (the prefix avoids
  masking base/stats functions such as `base::gamma` / `stats::poisson`): `dist_normal`
  (mean, **variance**), `dist_constant`, `dist_uniform`, `dist_gamma` (shape–scale),
  `dist_poisson`, `dist_beta`, `dist_exp`, `dist_lognormal`, and `dist_logistic`. They
  supply per-agent state-timer durations to the kernels, and `$sample_one()` / `$sample_n()`
  draw batches directly for interactive use and validation.

## Disease dynamics: the model menagerie

The kernels mutate the per-agent arrays and **return per-node counts**; the model applies
them to the census it maintains with `move_count(from, to, counts, tick)` (`from`/`to` may
be `NULL` for one-sided moves). So a model allocates only the states it has. All agent-loop
kernels parallelize across cores (Rayon) with private per-thread node accumulators.

- **`calc_foi(I, N, beta, seasonality, network, foi, tick)`** — per-node force of infection,
  with spatial coupling through a migration `network`. It reads the **settled
  start-of-interval** infectious census `I[t]`.
- **`transmission(...)`** — S → `to_state` (E or I), setting the `u16` timer;
  **`transmission_si(...)`** — S → I absorbing (no timer), the SI model.
- **Three step kernels** cover the whole menagerie, each a single `u16`-timer pass that
  branches on the agent's entry state and leads with M→S (so any model can add a maternal
  state): **`step_si`** (M→S, E→I), **`step_sir`** (+ I→absorbing S or R), **`step_sirs`**
  (+ I→R with an immunity timer, R→S).
- **Generic composable transitions** — `step_timer_expire(from, to)` and
  `step_timer_expire_set(from, to, duration)` express any single timed transition, the
  building blocks for models beyond the named menagerie (a vaccinated `V`, a second
  infectious stage, …) without writing Rust.

**One per-tick ordering for every model:** `carry_forward → step → calc_foi → transmission`,
with `calc_foi` immediately before `transmission` and no step kernel between them. Because
`calc_foi` reads the settled `I[t]`, an infectious agent contributes to the force of
infection on exactly the `D` census columns it occupies, so the effective reproduction
number is the **full `R0 = beta · D`** (never `beta · (D − 1)`) for both direct-S→I and
E-entry families, with no per-family special-casing. Validated against the Kermack–McKendrick
final-size relation `A = 1 − exp(−R0·A)` (SIR and SEIR attack fractions match to ~1e-3).

## Vital dynamics, demographics, and capacity

- **`births`** (newborns into a maternal `M` state, with a maternal timer and a
  Kaplan–Meier date of death), **`mortality`** (retire agents whose `dod` tick has arrived),
  **`constant_pop_vitals_sir`** (constant-population convenience), and **`import_infections`**
  (schedule-driven external seeding) grow and turn over an open population.
- **Realistic populations** — `AliasedDistribution` / `aliased_distribution` (Walker's alias
  method, O(1) per draw) sample ages from a pyramid; `KaplanMeierEstimator` /
  `kaplan_meier_estimator` draw a date of death conditioned on current age;
  `load_pyramid_csv` / `sample_pyramid_ages` wrap the laser.core pyramid format.
- **Sizing the agent array** — `calc_capacity` bounds the cumulative number ever born (use
  when slots are never reclaimed); `calc_capacity_cdr` bounds the **peak living** population
  for runs that reclaim dead slots with `squash`, so a century-long open run needs far fewer
  slots than one per agent ever born.

## Spatial coupling

- **`distances`** (haversine great-circle, km) plus the migration-network models
  **`gravity`**, **`radiation`**, **`stouffer`**, **`competing_destinations`**, and
  **`row_normalizer`** (caps each row's exported fraction) build the coupling matrix that
  `calc_foi` redistributes the force of infection through. Outputs match laser-core to
  floating-point precision.

## The `run_model()` runner

- **`run_model(scenario, model, nticks, r0, …)`** wires the whole eight-model menagerie
  (SI / SEI / SIS / SEIS / SIR / SEIR / SIRS / SEIRS) in the correct per-tick order, applies
  every `move_count` census delta, and returns a **`model` environment** bundling `$people`,
  `$nodes`, `$network`, `$carry`, `$states`, and the current `$tick`. It seeds any of the
  model's `E`/`I`/`R` states the scenario supplies and records the per-node flows
  `incidence`, `onset` (E→I), `recovery` (I-exit), and `waning` (R→S).
- **Extensibility** — optional spatial `network`, `seasonality`, and a reproducible `seed`;
  a `capacity` argument that reserves agent-array slots for `births` / `import_infections` to
  activate; an `extra_states` argument that registers states beyond S/E/I/R (a known name
  like `"M"` keeps its code and gets built-in M→S waning; a new name like `"V"` becomes a
  genuine new agent state the disease kernels leave untouched, driven entirely from
  callbacks); and the lifecycle callbacks `init(model)` plus per-tick `step_enter` /
  `step_update` (between the step kernel and `calc_foi`) / `step_exit`. These express
  constant-population vitals, importation, growth, vaccination, and quarantine without a
  hand-wired loop.
- **Validation & honesty** — rejects a non-finite/negative `r0` and `NA`/`Inf`/fractional/
  negative population and seed columns with clear messages, and warns when a scenario column
  or period argument is supplied that the chosen model does not use.

## Reproducibility, compaction, and utilities

- **`set_seed()` / `unset_seed()`** — after `set_seed(s)` an entire run is a deterministic
  function of `s` and the order of kernel calls, **independent of CPU/thread count** (the
  parallel kernels seed fixed-size agent chunks deterministically rather than per-thread).
  Without a seed, behaviour is entropy-seeded.
- **`squash(people)`** stably compacts living agents to the front of every per-agent Column
  (default mask: `state != D`) and updates the active count, reclaiming dead slots for reuse.
- **`values_map(value, n_ticks, n_nodes)`** broadcasts a scalar / per-node / per-tick /
  full-matrix value into a report Column (used to build `beta` / `seasonality` grids), and
  the **`bincount`** family — `bincount`, `bincount_wt` (weighted), `bincount_where`
  (predicate-filtered, count-aware), `bincount_where_wt` (weighted + filtered) — rolls
  per-agent properties up into per-node, per-tick reports in one parallel pass.

## Examples and documentation

- **Eight annotated teaching articles** on the [package website](https://clorton.github.io/razer/articles/)
  (`vignettes/articles/`, website-only) each pair 1:1 with a runnable script in `examples/`:
  *Getting started*, *Epidemic final size*, *Endemic dynamics*, *Spatial metapopulations*,
  *Demographics*, *Vital dynamics & measles*, *Interventions* (vaccination ± waning,
  quarantine), and *Long runs & memory*.
- Example scripts are **device-aware** (write to `examples/output/` under `Rscript`, draw to
  the Plots pane when `source()`d in RStudio); see `examples/README.md`.

## Infrastructure

- A **pkgdown website** (built in CI, deployed to `gh-pages`) with a curated function
  reference and the articles above. The changelog lives in `NEWS.md`.
- **GitHub Actions** for `R CMD check` and the pkgdown deploy, with every action pinned to a
  full commit SHA (version tag in a trailing comment) for supply-chain integrity.
