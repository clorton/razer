# Package index

## Population frame

Struct-of-arrays data store for agents and patches, and the epidemic
state codes used by the model kernels.

- [`LaserFrame`](https://clorton.github.io/razer/reference/LaserFrame.md)
  : Fixed-capacity struct-of-arrays population or patch data store.
- [`laser_states()`](https://clorton.github.io/razer/reference/laser_states.md)
  : Named integer vector of epidemic compartment state codes.

## Distributions

Parameterized probability distributions used to draw state-timer
durations (e.g. incubation and infectious periods).

- [`Distribution`](https://clorton.github.io/razer/reference/Distribution.md)
  : A parameterized probability distribution that can be sampled
  repeatedly.

- [`dist_beta()`](https://clorton.github.io/razer/reference/dist_beta.md)
  : Create a beta distribution on the open interval (0, 1).

- [`dist_constant()`](https://clorton.github.io/razer/reference/dist_constant.md)
  :

  Create a degenerate (constant) distribution that always returns
  `value`.

- [`dist_exp()`](https://clorton.github.io/razer/reference/dist_exp.md)
  :

  Create an exponential distribution with rate `rate`.

- [`dist_gamma()`](https://clorton.github.io/razer/reference/dist_gamma.md)
  : Create a gamma distribution parameterized by shape and scale.

- [`dist_logistic()`](https://clorton.github.io/razer/reference/dist_logistic.md)
  : Create a logistic distribution with the given location and scale.

- [`dist_lognormal()`](https://clorton.github.io/razer/reference/dist_lognormal.md)
  : Create a log-normal distribution.

- [`dist_normal()`](https://clorton.github.io/razer/reference/dist_normal.md)
  : Create a normal (Gaussian) distribution.

- [`dist_poisson()`](https://clorton.github.io/razer/reference/dist_poisson.md)
  :

  Create a Poisson distribution with rate (mean) `lambda`.

- [`dist_uniform()`](https://clorton.github.io/razer/reference/dist_uniform.md)
  : Create a continuous uniform distribution on the half-open interval
  \[low, high).

## Model step kernels

Per-tick transition kernels for composing SI / SIR / SEIR / SEIRS-family
models. Call them downstream-first within a tick — see “Composing
models: order of per-tick update operations” below.

- [`model-composition`](https://clorton.github.io/razer/reference/model-composition.md)
  [`update-order`](https://clorton.github.io/razer/reference/model-composition.md)
  : Composing models: order of per-tick update operations

- [`step_births_cbr()`](https://clorton.github.io/razer/reference/step_births_cbr.md)
  : Stochastic birth step using a crude birth rate.

- [`step_exposed_ei()`](https://clorton.github.io/razer/reference/step_exposed_ei.md)
  : Timer-based E→I transition (SEIR kernel).

- [`step_infectious_ir()`](https://clorton.github.io/razer/reference/step_infectious_ir.md)
  : Timer-based I→R transition (SIR / SEIR / SEIRS kernel).

- [`step_infectious_is()`](https://clorton.github.io/razer/reference/step_infectious_is.md)
  : Timer-based I→S transition for SIS models (no immunity).

- [`step_mortality_cdr()`](https://clorton.github.io/razer/reference/step_mortality_cdr.md)
  : Stochastic mortality step using a crude death rate.

- [`step_recovered_rs()`](https://clorton.github.io/razer/reference/step_recovered_rs.md)
  : Timer-based R→S transition for waning-immunity models.

- [`step_timer_expire()`](https://clorton.github.io/razer/reference/step_timer_expire.md)
  :

  Generalized timer-expiry transition into an *absorbing* (untimed)
  state.

- [`step_timer_expire_set()`](https://clorton.github.io/razer/reference/step_timer_expire_set.md)
  :

  Generalized timer-expiry transition into a state that has *its own*
  duration.

- [`step_transmission_se()`](https://clorton.github.io/razer/reference/step_transmission_se.md)
  : Stochastic S→E exposure step (SEIR kernel).

- [`step_transmission_si()`](https://clorton.github.io/razer/reference/step_transmission_si.md)
  : Stochastic S→I transmission step (SI / SIR kernel).
