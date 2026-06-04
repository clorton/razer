# razer examples

Runnable example scripts for the `razer` package. Install the package first
(`R CMD INSTALL .` or `devtools::install()` from the repository root), then run
an example with `Rscript`. Only base-R graphics are used (no extra dependencies).

Both examples run their transitions **downstream-first** within each tick (move
agents out of a compartment before moving agents in) so that exposed/infectious
periods are not artificially shortened by a tick — see the modeling note in the
repository `CLAUDE.md`.

## `sir_attack_fraction.R`

An agent-based SIR epidemic in a single well-mixed node, using
`step_infectious_ir` (I→R) then `step_transmission_si` (S→I) with a fixed
infectious period (`dist_constant`).

```sh
Rscript examples/sir_attack_fraction.R
```

Produces `sir_trajectories.png` (S/I/R over time) and `sir_attack_fraction.png`
(simulated final attack fraction across an `R0` sweep vs. the Kermack–McKendrick
final-size curve `A = 1 - exp(-R0 * A)`), and prints a comparison table plus
timing.

**Effective R0.** SIR uses `R0 = beta * (D - 1)`, where `D` is the infectious
period in ticks. An agent enters state `I` inside `step_transmission_si`, after
that step's infectious tally is computed, so it first contributes to the force of
infection on the *next* tick; it then recovers `D` ticks later. The net is `D - 1`
ticks of transmission. This one-tick offset is intrinsic to direct `S→I`
transmission and is independent of step order. With this mapping the simulated
attack fraction matches Kermack–McKendrick to within a few thousandths for
`R0 >= 1.5`.

## `seir_attack_fraction.R`

An agent-based SEIR epidemic, using `step_infectious_ir` (I→R), `step_exposed_ei`
(E→I), then `step_transmission_se` (S→E) with fixed incubation and infectious
periods.

```sh
Rscript examples/seir_attack_fraction.R
```

Produces `seir_trajectories.png` (S/E/I/R over time, showing the E peak preceding
the I peak) and `seir_attack_fraction.png` (the same Kermack–McKendrick
comparison).

**Effective R0.** SEIR uses `R0 = beta * D` (the *full* infectious period). Here an
agent enters `I` via `step_exposed_ei`, a separate step that runs *before* the
transmission tally in the downstream-first ordering, so it is counted in the force
of infection starting on its entry tick — recovering the full `D`-tick period. The
latent (exposed) period only delays the epidemic and does not affect the final
size, so the SEIR attack fraction obeys the same Kermack–McKendrick relation, and
the simulation matches it to within a few thousandths for `R0 >= 1.5`.

In both examples, just above the threshold (`R0` near 1) the finite seeded
population produces a small attack fraction where the deterministic theory predicts
zero — the expected near-critical behavior.
