# razer examples

Runnable example scripts for the `razer` package. Install the package first
(`R CMD INSTALL .` or `devtools::install()` from the repository root), then run
an example with `Rscript`.

## `sir_attack_fraction.R`

An agent-based SIR epidemic in a single well-mixed node, using
`step_transmission_si` (S→I) and `step_infectious_ir` (I→R) with a fixed
infectious period (`dist_constant`).

```sh
Rscript examples/sir_attack_fraction.R
```

It produces two plots in `examples/output/` and prints a comparison table:

- **`sir_trajectories.png`** — S / I / R counts over time for a single run.
- **`sir_attack_fraction.png`** — the simulated final attack fraction (fraction
  of the population ever infected, run to completion) across a sweep of the basic
  reproduction number `R0`, overlaid on the **Kermack–McKendrick** final-size
  curve `A = 1 - exp(-R0 * A)`.

Only base-R graphics are used (no extra package dependencies).

### Effective R0 in discrete time

The script maps the transmission rate to `R0 = beta * (D - 1)`, where `D` is the
infectious period in ticks. Each tick advances transmission first, then recovery:
a newly infected agent is given a timer of `D` during transmission and is then
decremented the same tick by `step_infectious_ir`, and it does not contribute to
that tick's force of infection. A secondary case therefore transmits on `D - 1`
ticks. With this mapping the simulated attack fraction matches the
Kermack–McKendrick prediction to within a few thousandths across the supercritical
range (`R0 >= 1.5`). Just above the threshold (`R0` near 1) the finite seeded
population produces a small attack fraction where the deterministic theory
predicts zero — the expected near-critical behavior.
