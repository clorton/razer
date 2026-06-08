# razer examples

Runnable example scripts for the `razer` package. Install the package first
(`R CMD INSTALL .` or `devtools::install()` from the repository root), then run an
example with `Rscript`. Only base-R graphics are used (no extra dependencies). Plots
are written to `examples/output/`.

All examples are built on the **Column kernels** (`allocate_scalar` / `allocate_vector`
buffers advanced by `calc_foi`, the `transmission` / `transmission_si` kernels, the
`step_si` / `step_sir` / `step_sirs` step kernels, `carry_forward_states` + `move_count`,
and the vital-dynamics kernels). The agent-loop kernels return per-node counts that the
model applies to its census. Every model uses one per-tick ordering (see the modeling note
in the repository `CLAUDE.md`): `carry_forward → step → calc_foi → transmission`, with
`calc_foi` immediately before `transmission`. `calc_foi` reads the settled start-of-interval
infectious census `I[t]`, so each infectious agent contributes to the force of infection on
exactly the `D` census columns it occupies — the effective reproduction number is the full
`R0 = beta * D` (never `beta * (D - 1)`), for both direct-S→I and E-entry models.

## `sir_attack_fraction.R`

An agent-based SIR epidemic in a single well-mixed node, advanced with the high-level
runner `run_model()` (which builds the population and runs `carry_forward_states` →
`step_sir` → `calc_foi` → `transmission` for you) with a fixed infectious period.

```sh
Rscript examples/sir_attack_fraction.R
```

Produces `sir_trajectories.png` (S/I/R over time) and `sir_attack_fraction.png`
(simulated final attack fraction across an `R0` sweep vs. the Kermack–McKendrick
final-size curve `A = 1 - exp(-R0 * A)`), and prints a comparison table plus timing.

**Effective R0 = `beta * D`.** `calc_foi` reads the settled start-of-interval infectious
census `I[t]`, so a directly-infected agent (which enters `I` one column after infection
and recovers `D` columns later) contributes to the force of infection on exactly its `D`
census columns. The simulated attack fraction matches Kermack–McKendrick to within a few
thousandths for `R0 >= 1.5`.

## `seir_attack_fraction.R`

An agent-based SEIR epidemic, advanced with `run_model()` (model `"SEIR"`), with fixed
incubation and infectious periods. The latent E stage only delays the epidemic; the final
size is unchanged.

```sh
Rscript examples/seir_attack_fraction.R
```

Produces `seir_trajectories.png` (S/E/I/R over time, showing the E peak preceding the I
peak) and `seir_attack_fraction.png` (the same Kermack–McKendrick comparison).

**Effective R0 = `beta * D`.** As in the SIR case, `calc_foi` reads the settled
start-of-interval census `I[t]`, so an agent (which enters `I` via the step kernel's E→I
one column after onset and recovers `D` columns later) contributes on exactly its `D`
census columns. The latent (exposed) period only delays the epidemic and does not affect
the final size, so the SEIR attack fraction obeys the same relation and matches it to
within a few thousandths for `R0 >= 1.5`.

In both attack-fraction examples, just above the threshold (`R0` near 1) the finite
seeded population produces a small attack fraction where the deterministic theory
predicts zero — the expected near-critical behavior.

## Other examples

All examples build their model through **`run_model()`**; the ones with vital dynamics,
importation, growth, or extra compartments add that behaviour through `run_model`'s
callbacks (`init` / `step_enter` / `step_update` / `step_exit`) and its `capacity` /
`extra_states` arguments, rather than a hand-wired per-tick loop.

- **`simple_sir.R`** — a spatial SIR over the England & Wales measles patches (radiation-
  model coupling network), run via `run_model("SIR", network=...)` with constant-population
  vital dynamics in a `step_exit` callback (`constant_pop_vitals_sir`).
- **`endemic_sir.R`** — a two-patch endemic SIR via `run_model()` with constant-pop vital
  turnover plus periodic importations (`step_exit` callback; `capacity` reserves slots for
  the imports); the susceptible fraction settles at `1/R0`.
- **`endemic_sir_seasonal.R`** — the endemic SIR with a gentle (±10%) annual sinusoid passed
  to `run_model`'s `seasonality`, producing phase-locked annual epidemic waves.
- **`engwal_measles.R`** — a full single-patch M-S-E-I-R measles model via
  `run_model("SEIR", extra_states = "M", capacity = ...)`: run_model tracks the maternal `M`
  compartment and applies its M→S waning, while births-into-M and Kaplan–Meier natural
  mortality run in a `step_exit` callback. Twenty years; recurrent epidemic waves.
- **`aged_population.R`** — builds a realistic age-structured population from a pyramid
  (alias sampler) and assigns each agent a Kaplan–Meier age at death.
- **`model_callbacks.R`** — extends `run_model()` through its `init` / `step_enter` /
  `step_exit` callbacks and the `model` environment: adds a per-agent date-of-birth Column,
  applies a one-time pulse vaccination mid-epidemic (S→R), and records a custom
  under-five-by-node report each tick with `bincount_where()`. Plots baseline vs.
  intervention.
- **`compare_models.R`** — runs `SIR`, `SIRS`, `SEIR`, and `SEIRS` through `run_model()` on
  the same population (1,000,000), duration (365 days), and transmission parameters
  (R0 = 2.5, infectious & incubation periods `gamma(140, 0.05)` ≈ 7 days), so only the
  compartment structure and waning differ. Trajectories are coloured by compartment
  (S blue, E orange, I red, R green) and styled by model (solid / dashed / dotted /
  long-dash); the waning models (SIRS, SEIRS) show recurrent epidemic waves while SIR/SEIR
  settle after one. Produces an overlay plot and a per-compartment 2×2 panel.
- **`long_run_squash.R`** — a 100-year demographic run (1,000,000 initial agents, CBR 30,
  CDR 15) that stays within a bounded agent array by reclaiming dead slots with `squash()`
  once a year. The array is sized with `calc_capacity_cdr()` (the peak-living bound, ~9.1M
  slots) instead of `calc_capacity()` (the cumulative-births bound, ~84M — 9× more). Plots
  the living population against both capacity bounds, the annual squash sawtooth in active
  slots, and daily births/deaths. Demonstrates `squash` + `calc_capacity_cdr` + the
  `births`/`mortality` vital-dynamics kernels end to end.
