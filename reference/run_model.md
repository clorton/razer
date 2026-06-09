# Run a closed-population agent-based model.

Builds an agent population from a per-node scenario, seeds it, and
advances one of the eight SI / SEI / SIS / SEIS / SIR / SEIR / SIRS /
SEIRS models for `nticks` ticks using the Column kernels in the correct
order (so the realized basic reproduction number is the full
`R0 = beta * D`). Each tick runs, uniformly for every model,
`carry_forward → step → calc_foi → transmission`, with `calc_foi` placed
immediately before `transmission` (no step kernel between them);
`run_model` applies every census delta for you.

## Usage

``` r
run_model(
  scenario,
  model,
  nticks,
  r0,
  infectious_period,
  incubation_period = NULL,
  immunity_period = NULL,
  network = NULL,
  seasonality = 1,
  seed = NULL,
  progress = FALSE,
  capacity = NULL,
  extra_states = NULL,
  init = NULL,
  step_enter = NULL,
  step_update = NULL,
  step_exit = NULL
)
```

## Arguments

- scenario:

  A data.frame with one row per node and an integer `population` column.
  Optional integer `E` / `I` / `R` columns seed that state per node —
  but only for states the model actually has (e.g. an `E` column is
  ignored by `SIR`). `S` is the remainder. The seeded `E + I + R` must
  not exceed a node's population.

- model:

  One of `"SI"`, `"SEI"`, `"SIS"`, `"SEIS"`, `"SIR"`, `"SEIR"`,
  `"SIRS"`, `"SEIRS"` (case-insensitive).

- nticks:

  Number of recorded daily states (dynamics run `nticks - 1` times).

- r0:

  Basic reproduction number; `beta = r0 / mean(infectious_period)`.

- infectious_period:

  Infectious period in ticks (a Distribution or a number).

- incubation_period:

  Latent/exposed period (required for the `SE*` models).

- immunity_period:

  Waning-immunity period (required for `SIRS` / `SEIRS`).

- network:

  Optional `n_nodes x n_nodes` spatial-coupling matrix (default: an
  all-zero matrix — independent nodes).

- seasonality:

  Transmission modifier in any
  [`values_map()`](https://clorton.github.io/razer/reference/values_map.md)-broadcastable
  form (default `1`).

- seed:

  Optional integer seed; if supplied,
  [`set_seed()`](https://clorton.github.io/razer/reference/set_seed.md)
  is called for a reproducible run.

- progress:

  If `TRUE`, draw a text progress bar.

- capacity:

  Optional agent-array capacity (\>= the total initial population, the
  default). Reserve extra slots — e.g.
  [`calc_capacity()`](https://clorton.github.io/razer/reference/calc_capacity.md)
  /
  [`calc_capacity_cdr()`](https://clorton.github.io/razer/reference/calc_capacity_cdr.md)
  — so a callback can grow the population with `births` or
  `import_infections` (which activate reserved slots). With the default
  `capacity`, the population is closed (no growth).

- extra_states:

  Optional character vector of extra agent-state names beyond the
  model's own S/E/I/R. A name already in
  [`laser_states()`](https://clorton.github.io/razer/reference/laser_states.md)
  (e.g. `"M"`) keeps that code; a NEW name (e.g. `"V"` for vaccinated)
  is assigned the next free code, becoming a genuine new agent state.
  Each gets a census Column that is carried forward, totalled into `N`,
  and seeded at tick 0 from the agents' states; the assigned codes are
  exposed in `model$states`. The disease kernels branch only on
  S/E/I/R/M, so an agent in a NEW state is left untouched (not infected,
  not transitioned) — its transitions are yours to drive from callbacks
  (move `S`→`V` to vaccinate; for waning, set a timer and run
  `step_timer_expire(V, S)` in `step_update`). For `"M"` only, run_model
  additionally applies the step kernels' built-in M→S waning each tick
  (recording the `waning_m` flow).

- init:

  Optional `function(model)` called once after the model is built and
  before the loop. Use it to add per-agent property
  [Column](https://clorton.github.io/razer/reference/Column.md)s, add a
  state census Column (append it to `model$carry` to have it carried
  forward), set agent states/timers, etc. The initial per-node census is
  (re)derived from the agents' states AFTER `init` runs, so changing
  agent states in `init` keeps the census consistent.

- step_enter:

  Optional `function(model)` called at the start of every tick (before
  the carry-forward). Receives the `model` environment; `model$tick` is
  the current 0-based interval index.

- step_update:

  Optional `function(model)` called every tick AFTER the step kernel
  (the disease-state progression) and BEFORE `calc_foi`. Because
  `calc_foi` reads the *settled* start-of-interval census `I[t]` (slice
  `model$tick`), edit the FOI **drivers** here — `model$nodes$beta`,
  `model$nodes$seasonality`, or `model$network` (all read at slice `t`)
  — to influence the CURRENT tick's force of infection. Census edits
  made here land in the working slice `t+1`, so they affect the NEXT
  tick's FOI, not this one.

- step_exit:

  Optional `function(model)` called at the end of every tick (after
  transmission). Receives the `model` environment (with `model$tick`
  set).

## Value

Invisibly, the `model` environment, with:

- `$people` — the agent arrays (`state`, `nodeid`, `timer`, plus
  `count`/`capacity`).

- `$nodes` — per-tick census Columns (`S`, `E`/`I`/`R` as applicable,
  any `extra_states` such as `M`), `N`, the `foi`, and the per-interval
  flows: `incidence` (new infections, all models), `onset` (E→I, `SE*`
  models), `recovery` (I to S or R, models with an I-exit), `waning`
  (R→S, `SIRS`/`SEIRS`), and `waning_m` (M→S, when `M` is an extra
  state) — each an `n_nodes`-wide, time-major Column.

- `$network` — the coupling weights as a 2-D f64
  [Column](https://clorton.github.io/razer/reference/Column.md)
  (`n_nodes x n_nodes`, column-major), built once from the `network`
  matrix and read by `calc_foi` each tick.

- `$carry` — the list of census Columns carried forward each tick.

- `$states` — the state-name→code map
  ([`laser_states()`](https://clorton.github.io/razer/reference/laser_states.md)
  plus any `extra_states`).

- `$tick` — during the loop, the current 0-based interval index (for the
  callbacks).

## Details

The `people`, `nodes`, `network`, and `carry` objects are bundled into a
single `model` environment, which is also what the optional callbacks
receive and what is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
scenario <- data.frame(population = 1e5L, I = 100L)
m <- run_model(scenario, "SEIR", nticks = 200L, r0 = 2.5,
               infectious_period = 8, incubation_period = 5, seed = 1)
tail(m$nodes$I$values())
} # }
```
