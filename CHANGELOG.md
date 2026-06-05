# Changelog

All notable changes to this project are documented here.

## Unreleased

### Changed
- **Breaking:** `step_transmission_si()` and `step_transmission_se()` now take a
  required fifth argument, `network`, enabling spatial coupling of contagion (see
  Added). razer is spatial-first: there is no "no coupling" special case — pass an
  **all-zero matrix** to run nodes independently (identical to the previous
  per-node force of infection). All existing callers were updated.

### Added
- `bincount(values, nbins, counts)` — a parallel, NumPy-`bincount`-style histogram
  over a [Column] (`src/rust/src/bincount.rs`). For each bin `b` in `0..nbins` it
  counts how many elements of the integer-typed `values` Column equal `b` and
  writes the result into the caller-provided `counts` Column (length `>= nbins`;
  entries at/beyond `nbins` are left untouched). Each Rayon worker accumulates into
  a private per-thread histogram (no shared-bin write collisions); the used range
  of `counts` is then zeroed and the per-thread tallies are reduced into it. One
  generic kernel serves every value width/signedness (`i8`..`u32`) and any numeric
  `counts` type via traits. Counts 30.18M agents into 954 node bins in ~15 ms;
  covered by `tests/testthat/test-bincount.R` (incl. a parallel-vs-`tabulate`
  cross-check).
- `bincountw(values, weights, nbins, counts)` — the weighted counterpart: sums
  each element's weight into its bin (`counts[b] = Σ weights[i]` over `i` with
  `values[i] == b`), à la `numpy.bincount(values, weights=...)`. Same parallel,
  collision-free, zero-then-accumulate design; `weights` may be any numeric
  `Column` (signed, unsigned, or floating point — all widened to f64 for
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
  before transmission S→I). It records per-node compartment trajectories (`S`,
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
  runnable SIR and SEIR examples that plot the compartment trajectories and
  compare the simulated final attack fraction against the Kermack–McKendrick
  final-size relation `A = 1 - exp(-R0 * A)` across an `R0` sweep, with timing
  output (with `examples/README.md` and sample output plots).
- `CLAUDE.md` documenting the downstream-first transition-ordering convention for
  composing models from the step kernels.

### Changed
- Model compositions (SIR/SEIR/SEIRS test helpers and the examples) now call the
  per-tick transitions **downstream-first** (out of each timed compartment before
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
