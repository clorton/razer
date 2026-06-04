#!/usr/bin/env Rscript

# SIR example: trajectories and the Kermack–McKendrick final-size relation.
#
# This script
#   1. runs razer's agent-based SIR kernel for a single parameter set and plots
#      the S / I / R trajectories over time, then
#   2. sweeps the basic reproduction number R0 and compares the *simulated*
#      attack fraction (the fraction of the population ever infected) against the
#      classic Kermack–McKendrick final-size equation
#
#          A = 1 - exp(-R0 * A)
#
#      whose non-trivial root gives the attack fraction for a large, well-mixed,
#      fully susceptible population.
#
# Transition ordering: each tick runs transitions downstream-first (recovery
# I->R before transmission S->I), so a newly infected agent is never decremented
# in the same tick it is infected. See the package modeling note for why this
# ordering matters in general.
#
# Discrete-time mapping: with per-tick force of infection p = 1 - exp(-beta*I/N),
# an infectious agent infects on the order of beta susceptibles per tick while the
# population is still mostly susceptible. The subtlety is *how many* ticks a
# secondary case transmits for. In SIR the agent enters state I inside
# step_transmission_si, which computes its infectious tally at the start of the
# step — before that step's new infections — so a newly infected agent first
# contributes to the force of infection on the *next* tick. It then recovers after
# D ticks, and on its recovery tick it is moved to R (by step_infectious_ir, which
# runs first) before the tally is taken. The net effect is that a secondary case
# is counted in the force of infection on D - 1 ticks, giving an effective basic
# reproduction number
#
#     R0 = beta * (D - 1).
#
# (This one-tick offset is specific to direct S->I transmission, where infection
# and the infectious tally share a step. In SEIR an agent enters I via a separate
# step, so the full period D applies — see examples/seir_attack_fraction.R.)
#
# With this mapping the deterministic final size satisfies the equation above
# exactly, provided the epidemic is run to completion — which the sweep below does
# by iterating until no infectious agents remain.
#
# Run from anywhere:  Rscript examples/sir_attack_fraction.R
# Output PNGs are written next to this script in examples/output/.

library(razer)

set.seed(20260604L)

states <- laser_states()   # c(S = 0, E = 1, I = 2, R = 3, D = -1)

# ── build + seed a single-node SIR population ───────────────────────────────────
new_sir <- function(n, n_seed, inf_duration) {
  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", states[["S"]])
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)

  # Seed the first n_seed agents as infectious with a full infectious timer.
  state0 <- rep(states[["S"]], n)
  state0[seq_len(n_seed)] <- states[["I"]]
  ppl$state <- state0
  timer0 <- rep(0L, n)
  timer0[seq_len(n_seed)] <- inf_duration
  ppl$timer <- timer0

  # Single well-mixed node holding the total population and an I tally.
  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  list(ppl = ppl, nd = nd,
       inf_dist = dist_constant(inf_duration),  # fixed infectious period
       imm_dist = dist_constant(0))             # SIR: no waning, R timer unused
}

# Advance one tick. Transitions run downstream-first: I->R recovery before S->I
# transmission, so an agent infected this tick is not decremented by the recovery
# step in the same tick (see header note on effective infectious period).
sir_step <- function(sim, beta) {
  step_infectious_ir(sim$ppl, imm_dist = sim$imm_dist)
  step_transmission_si(sim$ppl, sim$nd, beta = beta, inf_dist = sim$inf_dist)
}

n_in_state <- function(sim, code) sum(sim$ppl$state == code)

# ── one SIR run over a fixed horizon, returning the S / I / R trajectory ─────────
run_sir_traj <- function(n, n_seed, beta, inf_duration, nticks) {
  sim <- new_sir(n, n_seed, inf_duration)
  traj <- matrix(0L, nrow = nticks + 1L, ncol = 3L,
                 dimnames = list(NULL, c("S", "I", "R")))
  traj[1L, ] <- c(n - n_seed, n_seed, 0L)
  for (tick in seq_len(nticks)) {
    sir_step(sim, beta)
    traj[tick + 1L, ] <- c(n_in_state(sim, states[["S"]]),
                           n_in_state(sim, states[["I"]]),
                           n_in_state(sim, states[["R"]]))
  }
  traj
}

# ── run to completion (no infectious agents left) and return the attack fraction ─
final_attack_fraction <- function(n, n_seed, beta, inf_duration, max_ticks = 5000L) {
  sim <- new_sir(n, n_seed, inf_duration)
  for (tick in seq_len(max_ticks)) {
    sir_step(sim, beta)
    if (n_in_state(sim, states[["I"]]) == 0L) break
  }
  (n - n_in_state(sim, states[["S"]])) / n   # fraction ever infected
}

# ── Kermack–McKendrick final size ──────────────────────────────────────────────
# Solve A = 1 - exp(-R0 * A) for the non-trivial root in (0, 1). For R0 <= 1 the
# only root is A = 0 (no epidemic).
km_attack_fraction <- function(R0) {
  if (R0 <= 1) return(0)
  f <- function(A) 1 - exp(-R0 * A) - A
  uniroot(f, lower = 1e-9, upper = 1 - 1e-12)$root
}

# ── output location (next to this script) ───────────────────────────────────────
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N <- 100000L
D <- 10L

# ── 1. single-run trajectories ──────────────────────────────────────────────────
beta_traj <- 0.25
R0_traj   <- beta_traj * (D - 1)        # effective R0 (see header)
nticks_traj <- 200L
traj_time <- system.time(
  traj <- run_sir_traj(N, n_seed = 100L, beta = beta_traj, inf_duration = D, nticks = nticks_traj)
)
cat(sprintf("run_sir_traj: %s agents x %d ticks took %.3f s (%.1f ns/agent-tick)\n",
            format(N, big.mark = ","), nticks_traj, traj_time[["elapsed"]],
            1e9 * traj_time[["elapsed"]] / (as.numeric(N) * nticks_traj)))

png(file.path(out_dir, "sir_trajectories.png"), width = 1000, height = 650, res = 110)
ticks <- 0:(nrow(traj) - 1L)
matplot(ticks, traj, type = "l", lwd = 2.5, lty = 1,
        col = c("#2c7fb8", "#d7301f", "#238b45"),
        xlab = "tick", ylab = "number of agents",
        main = sprintf("razer SIR trajectories  (N = %s, R0 = beta*(D-1) = %.2f)",
                       format(N, big.mark = ","), R0_traj))
legend("right", legend = c("S", "I", "R"), bty = "n", lwd = 2.5,
       col = c("#2c7fb8", "#d7301f", "#238b45"))
dev.off()

# ── 2. attack-fraction sweep vs Kermack–McKendrick ───────────────────────────────
R0_grid <- seq(0.5, 4.0, by = 0.25)
sweep_time <- system.time(
  sim_af <- vapply(R0_grid, function(R0)
    final_attack_fraction(N, n_seed = 50L, beta = R0 / (D - 1), inf_duration = D),
    numeric(1))
)
cat(sprintf("attack-fraction sweep: %d runs to completion (N = %s each) took %.3f s\n",
            length(R0_grid), format(N, big.mark = ","), sweep_time[["elapsed"]]))
km_af   <- vapply(R0_grid, km_attack_fraction, numeric(1))

png(file.path(out_dir, "sir_attack_fraction.png"), width = 1000, height = 650, res = 110)
plot(R0_grid, km_af, type = "l", lwd = 2.5, col = "black",
     ylim = c(0, 1), xlab = expression(R[0]), ylab = "attack fraction",
     main = "Attack fraction: razer SIR vs Kermack–McKendrick")
points(R0_grid, sim_af, pch = 19, cex = 1.1, col = "#d7301f")
abline(v = 1, lty = 3, col = "grey50")
legend("bottomright", bty = "n",
       legend = c("Kermack-McKendrick  A = 1 - exp(-R0 A)", "razer simulation"),
       lwd = c(2.5, NA), pch = c(NA, 19), col = c("black", "#d7301f"))
dev.off()

# ── comparison table ─────────────────────────────────────────────────────────────
cmp <- data.frame(
  R0                 = R0_grid,
  simulated          = round(sim_af, 4),
  kermack_mckendrick = round(km_af, 4),
  abs_error          = round(abs(sim_af - km_af), 4)
)
cat("\nAttack fraction: simulated vs Kermack-McKendrick\n")
print(cmp, row.names = FALSE)
cat(sprintf("\nMax absolute error for R0 >= 1.5: %.4f\n",
            max(cmp$abs_error[cmp$R0 >= 1.5])))
cat(sprintf("Plots written to: %s\n", normalizePath(out_dir)))
