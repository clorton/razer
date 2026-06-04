# Composing models: order of per-tick update operations

Guidance for assembling a model from the per-tick step kernels
([`step_transmission_si()`](https://clorton.github.io/razer/reference/step_transmission_si.md),
[`step_transmission_se()`](https://clorton.github.io/razer/reference/step_transmission_se.md),
[`step_exposed_ei()`](https://clorton.github.io/razer/reference/step_exposed_ei.md),
[`step_infectious_ir()`](https://clorton.github.io/razer/reference/step_infectious_ir.md),
[`step_infectious_is()`](https://clorton.github.io/razer/reference/step_infectious_is.md),
[`step_recovered_rs()`](https://clorton.github.io/razer/reference/step_recovered_rs.md)):
within each tick, call the transitions **downstream-first** — move
agents *out* of a timed compartment before moving agents *in*, walking
the disease-progression chain backwards.

## Details

Recommended per-tick call order, by model:

- SIR:

  [`step_infectious_ir()`](https://clorton.github.io/razer/reference/step_infectious_ir.md),
  then
  [`step_transmission_si()`](https://clorton.github.io/razer/reference/step_transmission_si.md)

- SEIR:

  [`step_infectious_ir()`](https://clorton.github.io/razer/reference/step_infectious_ir.md),
  [`step_exposed_ei()`](https://clorton.github.io/razer/reference/step_exposed_ei.md),
  [`step_transmission_se()`](https://clorton.github.io/razer/reference/step_transmission_se.md)

- SEIRS:

  [`step_recovered_rs()`](https://clorton.github.io/razer/reference/step_recovered_rs.md),
  [`step_infectious_ir()`](https://clorton.github.io/razer/reference/step_infectious_ir.md),
  [`step_exposed_ei()`](https://clorton.github.io/razer/reference/step_exposed_ei.md),
  [`step_transmission_se()`](https://clorton.github.io/razer/reference/step_transmission_se.md)

- SIS:

  [`step_infectious_is()`](https://clorton.github.io/razer/reference/step_infectious_is.md),
  then
  [`step_transmission_si()`](https://clorton.github.io/razer/reference/step_transmission_si.md)

**Why it matters.** Each timed transition sets a fresh duration timer
when an agent *enters* a state, and the step for the *next* transition
decrements that timer. If the inflow step runs before the outflow step
within the same tick, a just-arrived agent is decremented immediately
and loses one tick of residency — its effective exposed or infectious
period becomes `duration - 1` rather than `duration`. Running
downstream-first prevents this: a newly exposed or infectious agent is
not touched again until the next tick.

**SIR caveat.** For direct S -\> I transmission a one-tick offset
persists regardless of order. An agent becomes infectious *inside*
[`step_transmission_si()`](https://clorton.github.io/razer/reference/step_transmission_si.md),
after that step has already computed its infectious tally, so it first
contributes to the force of infection on the *next* tick. The SIR
effective basic reproduction number is therefore `beta * (D - 1)`,
whereas SEIR — where agents enter `I` via the separate
[`step_exposed_ei()`](https://clorton.github.io/razer/reference/step_exposed_ei.md),
which runs before the transmission tally — yields the full `beta * D`.

The `sir_attack_fraction.R` and `seir_attack_fraction.R` scripts in the
package's `examples/` directory demonstrate both, validating the
simulated attack fraction against the Kermack-McKendrick final-size
relation.

## See also

[`step_transmission_si()`](https://clorton.github.io/razer/reference/step_transmission_si.md),
[`step_transmission_se()`](https://clorton.github.io/razer/reference/step_transmission_se.md),
[`step_exposed_ei()`](https://clorton.github.io/razer/reference/step_exposed_ei.md),
[`step_infectious_ir()`](https://clorton.github.io/razer/reference/step_infectious_ir.md),
[`step_infectious_is()`](https://clorton.github.io/razer/reference/step_infectious_is.md),
[`step_recovered_rs()`](https://clorton.github.io/razer/reference/step_recovered_rs.md)
