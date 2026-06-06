# razer examples

Runnable example scripts for the `razer` package. Install the package first
(`R CMD INSTALL .` or `devtools::install()` from the repository root), then run an
example with `Rscript`. Only base-R graphics are used (no extra dependencies). Plots
are written to `examples/output/`.

All examples are built on the **Column kernels** (`allocate_scalar` / `allocate_vector`
buffers advanced by `calc_foi`, `sir_step`, `measles_step`, `transmission` /
`transmission_u16`, `carry_forward_states`, and the vital-dynamics kernels). The
per-tick ordering follows the modeling note in the repository `CLAUDE.md`; in
particular `calc_foi` is placed so the effective reproduction number is the full
`R0 = beta * D` (never `beta * (D - 1)`).

## `sir_attack_fraction.R`

An agent-based SIR epidemic in a single well-mixed node. Per tick:
`carry_forward_states` → `calc_foi` → `sir_step` (I→R) → `transmission` (S→I) with a
fixed infectious period (`dist_constant`).

```sh
Rscript examples/sir_attack_fraction.R
```

Produces `sir_trajectories.png` (S/I/R over time) and `sir_attack_fraction.png`
(simulated final attack fraction across an `R0` sweep vs. the Kermack–McKendrick
final-size curve `A = 1 - exp(-R0 * A)`), and prints a comparison table plus timing.

**Effective R0 = `beta * D`.** A directly-infected agent enters `I` *after* the tick's
force-of-infection tally, so it is not counted on its entry tick. Running `calc_foi`
*before* the recovery step (`sir_step`) counts it on its recovery tick instead — losing
the entry tick but gaining the recovery tick nets the full `D`-tick infectious period.
The simulated attack fraction matches Kermack–McKendrick to within a few thousandths for
`R0 >= 1.5`.

## `seir_attack_fraction.R`

An agent-based SEIR epidemic (the measles kernels with the maternal compartment M left
empty). Per tick: `carry_forward_states` → `measles_step` (E→I, I→R) → `calc_foi` →
`transmission_u16` (S→E) with fixed incubation and infectious periods.

```sh
Rscript examples/seir_attack_fraction.R
```

Produces `seir_trajectories.png` (S/E/I/R over time, showing the E peak preceding the I
peak) and `seir_attack_fraction.png` (the same Kermack–McKendrick comparison).

**Effective R0 = `beta * D`.** An agent enters `I` via `measles_step`, a step run
*before* the force-of-infection tally, so it is counted from its entry tick — the full
`D`-tick period. The latent (exposed) period only delays the epidemic and does not affect
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
