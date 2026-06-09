# Package index

## Model runner

One call to build, seed, and run a closed-population agent-based model
(SI / SEI / SIS / SEIS / SIR / SEIR / SIRS / SEIRS) on the Column
kernels.

- [`run_model()`](https://clorton.github.io/razer/reference/run_model.md)
  : Run a closed-population agent-based model.

## Agent & node storage

Rust-owned, dtype-tagged property arrays (per-agent and
per-node/per-tick), the incremental-census carry-forward helpers, and
the epidemic state codes.

- [`Column`](https://clorton.github.io/razer/reference/Column.md) : A
  Rust-owned, dtype-tagged property array (1-D scalar, or a 2-D vector
  report).

- [`allocate_scalar()`](https://clorton.github.io/razer/reference/allocate_scalar.md)
  : Allocate a fresh, zero-filled property array of a given type and
  length.

- [`allocate_vector()`](https://clorton.github.io/razer/reference/allocate_vector.md)
  : Allocate a fresh, zero-filled 2-D property array (a per-slot report
  buffer).

- [`carry_forward()`](https://clorton.github.io/razer/reference/carry_forward.md)
  :

  Carry a per-node counter forward one tick: copy column `tick` onto
  `tick + 1`.

- [`carry_forward_states()`](https://clorton.github.io/razer/reference/carry_forward_states.md)
  : Carry census counters forward, and optionally total some of them.

- [`laser_states()`](https://clorton.github.io/razer/reference/laser_states.md)
  : Named integer vector of epidemic state codes.

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

## Transmission & disease dynamics

Per-tick Column kernels for force of infection and the timed state
transitions of the SI/SEI/SIS/SEIS/SIR/SEIR/SIRS/SEIRS menagerie. The
kernels return per-node counts the model applies to its census via
`move_count`; place `calc_foi` so the effective reproduction number is
the full `beta * D` (see the modeling note in the package README /
CLAUDE.md).

- [`calc_foi()`](https://clorton.github.io/razer/reference/calc_foi.md)
  : Compute the per-node force of infection (FOI) for one tick.

- [`transmission()`](https://clorton.github.io/razer/reference/transmission.md)
  :

  Stochastic transmission S→`to_state`, returning new infections per
  node.

- [`transmission_si()`](https://clorton.github.io/razer/reference/transmission_si.md)
  :

  Stochastic transmission S→I into an ABSORBING `I` (the SI model),
  returning new infections per node.

- [`step_si()`](https://clorton.github.io/razer/reference/step_si.md) :
  Advance M→S (maternal waning) and E→I (incubation) for one tick — SI /
  SEI.

- [`step_sir()`](https://clorton.github.io/razer/reference/step_sir.md)
  :

  Advance M→S, E→I, and I→`absorbing_state` for one tick — SIS / SIR /
  SEIS / SEIR.

- [`step_sirs()`](https://clorton.github.io/razer/reference/step_sirs.md)
  : Advance M→S, E→I, I→R (with waning immunity), and R→S for one tick —
  SIRS / SEIRS.

- [`step_timer_expire()`](https://clorton.github.io/razer/reference/step_timer_expire.md)
  :

  Generic timed transition `from_state -> to_state` into an UNTIMED
  destination.

- [`step_timer_expire_set()`](https://clorton.github.io/razer/reference/step_timer_expire_set.md)
  :

  Generic timed transition `from_state -> to_state` into a TIMED
  destination.

- [`move_count()`](https://clorton.github.io/razer/reference/move_count.md)
  :

  Apply a per-node transition count to the census at `tick + 1`.

## Vital dynamics & importation

Births, deaths, and externally-seeded cases.

- [`births()`](https://clorton.github.io/razer/reference/births.md) :
  Apply crude-birth-rate births for one tick; newborns enter maternal
  immunity (M).
- [`mortality()`](https://clorton.github.io/razer/reference/mortality.md)
  : Apply natural mortality for one tick, returning deaths per node by
  state.
- [`constant_pop_vitals_sir()`](https://clorton.github.io/razer/reference/constant_pop_vitals_sir.md)
  : Apply constant-population SIR vital dynamics for one tick.
- [`import_infections()`](https://clorton.github.io/razer/reference/import_infections.md)
  : Import new infectious cases from a schedule, activating reserved
  agent slots.

## Demographics & initialization

Realistic age structure and dates of death, and capacity sizing.

- [`AliasedDistribution`](https://clorton.github.io/razer/reference/AliasedDistribution.md)
  :

  A discrete distribution over bin indices `0..n`, sampled by the Vose
  alias method.

- [`aliased_distribution()`](https://clorton.github.io/razer/reference/aliased_distribution.md)
  : Build an AliasedDistribution from a vector of non-negative bin
  counts.

- [`KaplanMeierEstimator`](https://clorton.github.io/razer/reference/KaplanMeierEstimator.md)
  : A Kaplan–Meier sampler over a cumulative-deaths-by-year life table.

- [`kaplan_meier_estimator()`](https://clorton.github.io/razer/reference/kaplan_meier_estimator.md)
  : Build a KaplanMeierEstimator from cumulative deaths by year.

- [`load_pyramid_csv()`](https://clorton.github.io/razer/reference/load_pyramid_csv.md)
  : Load a population-pyramid CSV into a numeric matrix.

- [`sample_pyramid_ages()`](https://clorton.github.io/razer/reference/sample_pyramid_ages.md)
  : Sample realistic per-agent ages (in days) from a population pyramid.

- [`calc_capacity()`](https://clorton.github.io/razer/reference/calc_capacity.md)
  : Estimate the agent capacity needed for a growing population.

- [`calc_capacity_cdr()`](https://clorton.github.io/razer/reference/calc_capacity_cdr.md)
  :

  Estimate the agent capacity for a growing population reclaimed with
  [`squash()`](https://clorton.github.io/razer/reference/squash.md).

## Spatial coupling & migration

Pairwise distances and migration-network models for spatial coupling of
the force of infection.

- [`distances()`](https://clorton.github.io/razer/reference/distances.md)
  : Great-circle distance matrix between geographic points.

- [`gravity()`](https://clorton.github.io/razer/reference/gravity.md) :
  Gravity migration-network model.

- [`radiation()`](https://clorton.github.io/razer/reference/radiation.md)
  : Radiation migration-network model (Simini et al., Nature 2012).

- [`stouffer()`](https://clorton.github.io/razer/reference/stouffer.md)
  : Stouffer's intervening-opportunities migration model (Stouffer,
  1940).

- [`competing_destinations()`](https://clorton.github.io/razer/reference/competing_destinations.md)
  : Competing-destinations migration model (Fotheringham, 1984).

- [`row_normalizer()`](https://clorton.github.io/razer/reference/row_normalizer.md)
  :

  Cap each row sum of a network matrix at `max_rowsum`.

## Reproducibility

Seed the simulation RNG. After
[`set_seed()`](https://clorton.github.io/razer/reference/set_seed.md), a
run is reproducible regardless of CPU/thread count;
[`unset_seed()`](https://clorton.github.io/razer/reference/unset_seed.md)
reverts to entropy-seeded (random) runs.

- [`set_seed()`](https://clorton.github.io/razer/reference/set_seed.md)
  : Set the global random seed, making subsequent razer runs
  reproducible.

- [`unset_seed()`](https://clorton.github.io/razer/reference/unset_seed.md)
  :

  Revert to a non-reproducible, entropy-seeded RNG (undo
  [`set_seed()`](https://clorton.github.io/razer/reference/set_seed.md)).

## Utilities

Grid broadcasting and binning helpers.

- [`values_map()`](https://clorton.github.io/razer/reference/values_map.md)
  :

  Build a values map: broadcast a value into a `n_ticks x n_nodes` grid
  Column.

- [`bincount()`](https://clorton.github.io/razer/reference/bincount.md)
  :

  Count occurrences of each value, NumPy `bincount`-style, into a
  buffer.

- [`bincount_wt()`](https://clorton.github.io/razer/reference/bincount_wt.md)
  : Weighted bincount: sum each element's weight into its bin.

- [`bincount_where()`](https://clorton.github.io/razer/reference/bincount_where.md)
  : Count, per group, the agents whose property satisfies a comparison.

- [`bincount_where_wt()`](https://clorton.github.io/razer/reference/bincount_where_wt.md)
  : Sum a weight, per group, over the agents whose property satisfies a
  comparison.

- [`squash()`](https://clorton.github.io/razer/reference/squash.md) :
  Compact a people environment, dropping excluded agents and reclaiming
  their slots.
