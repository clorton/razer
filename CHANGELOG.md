# Changelog

All notable changes to this project are documented here.

## Unreleased

### Added
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
- `examples/sir_attack_fraction.R` — a runnable SIR example that plots the S/I/R
  trajectories and compares the simulated final attack fraction against the
  Kermack–McKendrick final-size relation `A = 1 - exp(-R0 * A)` across an `R0`
  sweep (with `examples/README.md` and sample output plots).

### Changed
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
