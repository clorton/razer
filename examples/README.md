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

An agent-based SIR epidemic in a single well-mixed node. Per tick:
`carry_forward_states` → `step_sir` (absorbing = R; I→R) → `calc_foi` → `transmission`
(S→I) with a fixed infectious period (`dist_constant`).

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

An agent-based SEIR epidemic (`step_sir` with absorbing = R; no M compartment). Per tick:
`carry_forward_states` → `step_sir` (E→I, I→R) → `calc_foi` → `transmission` (S→E) with
fixed incubation and infectious periods.

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

- **`simple_sir.R`** — a spatial SIR over the England & Wales measles patches: builds a
  distance matrix, a radiation-model coupling network, seeds infections, and runs the
  Column SIR loop with constant-population vital dynamics.
- **`endemic_sir.R`** — a two-patch endemic SIR sustained by vital turnover plus periodic
  importations; the susceptible fraction settles at `1/R0`.
- **`endemic_sir_seasonal.R`** — the endemic SIR with a gentle (±10%) annual sinusoid on
  transmission, producing phase-locked annual epidemic waves.
- **`engwal_measles.R`** — a full single-patch M-S-E-I-R measles model with maternal
  immunity, a Kaplan–Meier date of death per agent, natural mortality, births into M, and
  a CBR-sized agent capacity, run for 20 years.
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
