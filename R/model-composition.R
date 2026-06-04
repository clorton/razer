#' Composing models: order of per-tick update operations
#'
#' @description
#' Guidance for assembling a model from the per-tick step kernels
#' ([step_transmission_si()], [step_transmission_se()], [step_exposed_ei()],
#' [step_infectious_ir()], [step_infectious_is()], [step_recovered_rs()]): within
#' each tick, call the transitions **downstream-first** — move agents *out* of a
#' timed compartment before moving agents *in*, walking the disease-progression
#' chain backwards.
#'
#' @details
#' Recommended per-tick call order, by model:
#'
#' \describe{
#'   \item{SIR}{`step_infectious_ir()`, then `step_transmission_si()`}
#'   \item{SEIR}{`step_infectious_ir()`, `step_exposed_ei()`, `step_transmission_se()`}
#'   \item{SEIRS}{`step_recovered_rs()`, `step_infectious_ir()`, `step_exposed_ei()`, `step_transmission_se()`}
#'   \item{SIS}{`step_infectious_is()`, then `step_transmission_si()`}
#' }
#'
#' **Why it matters.** Each timed transition sets a fresh duration timer when an
#' agent *enters* a state, and the step for the *next* transition decrements that
#' timer. If the inflow step runs before the outflow step within the same tick, a
#' just-arrived agent is decremented immediately and loses one tick of residency —
#' its effective exposed or infectious period becomes `duration - 1` rather than
#' `duration`. Running downstream-first prevents this: a newly exposed or
#' infectious agent is not touched again until the next tick.
#'
#' **SIR caveat.** For direct S -> I transmission a one-tick offset persists
#' regardless of order. An agent becomes infectious *inside*
#' [step_transmission_si()], after that step has already computed its infectious
#' tally, so it first contributes to the force of infection on the *next* tick.
#' The SIR effective basic reproduction number is therefore `beta * (D - 1)`,
#' whereas SEIR — where agents enter `I` via the separate [step_exposed_ei()],
#' which runs before the transmission tally — yields the full `beta * D`.
#'
#' The `sir_attack_fraction.R` and `seir_attack_fraction.R` scripts in the
#' package's `examples/` directory demonstrate both, validating the simulated
#' attack fraction against the Kermack-McKendrick final-size relation.
#'
#' **Generalized timer kernels.** The four named timer transitions are
#' fixed-state shorthands over two generalized kernels parameterized by the
#' `from`/`to` state codes — use these directly to add transitions the named
#' kernels do not cover:
#' \describe{
#'   \item{[step_timer_expire()]`(people, from, to)`}{transition into an
#'     *absorbing* (untimed) state; the timer is left at 0. Generalizes
#'     [step_infectious_is()] (I->S) and [step_recovered_rs()] (R->S).}
#'   \item{[step_timer_expire_set()]`(people, from, to, dist)`}{transition into a
#'     state with *its own* duration; a fresh per-agent timer is drawn from
#'     `dist`. Generalizes [step_exposed_ei()] (E->I) and [step_infectious_ir()]
#'     (I->R).}
#' }
#'
#' @seealso [step_transmission_si()], [step_transmission_se()],
#'   [step_timer_expire()], [step_timer_expire_set()],
#'   [step_exposed_ei()], [step_infectious_ir()], [step_infectious_is()],
#'   [step_recovered_rs()]
#' @name model-composition
#' @aliases model-composition update-order
#' @keywords documentation
NULL
