#!/usr/bin/env Rscript

# SEIR example: trajectories and the Kermack–McKendrick final-size relation.
#
# Companion to sir_attack_fraction.R. It
#   1. runs razer's agent-based SEIR kernel for a single parameter set and plots
#      the S / E / I / R trajectories over time, then
#   2. sweeps the basic reproduction number R0 and compares the simulated attack
#      fraction (fraction of the population ever infected) against the
#      Kermack–McKendrick final-size equation
#
#          A = 1 - exp(-R0 * A).
#
# The latent (exposed) period only *delays* the epidemic; it does not change the
# final size, so the SEIR attack fraction obeys the same relation as SIR with
# R0 determined by the infectious period alone.
#
# Transition ordering: each tick runs transitions downstream-first --
# I->R, then E->I, then S->E -- so an agent entering E or I is not decremented in
# the same tick it arrives. Because an agent enters I via step_exposed_ei (a
# separate step that runs *before* the transmission tally), it is counted in the
# force of infection starting on its entry tick, so the effective infectious
# period is the full D ticks and
#
#     R0 = beta * D.
#
# (Contrast SIR, where infection and the infectious tally share one step, giving
# R0 = beta * (D - 1); see sir_attack_fraction.R.)
#
# Run from anywhere:  Rscript examples/seir_attack_fraction.R
# Output PNGs are written next to this script in examples/output/.

library(razer)

set.seed(20260604L)

states <- laser_states()   # c(S = 0, E = 1, I = 2, R = 3, D = -1)

# ── build + seed a single-node SEIR population ──────────────────────────────────
new_seir <- function(n, n_seed, exp_duration, inf_duration) {
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

  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  list(ppl = ppl, nd = nd,
       exp_dist = dist_constant(exp_duration),  # incubation period
       inf_dist = dist_constant(inf_duration),  # infectious period
       imm_dist = dist_constant(0))             # SEIR: no waning, R timer unused
}

# Advance one tick, downstream-first: I->R, then E->I, then S->E. This keeps the
# exposed and infectious periods at their full durations (a newly arrived E or I
# agent is not decremented by the next step in the same tick).
seir_step <- function(sim, beta) {
  step_infectious_ir(sim$ppl, imm_dist = sim$imm_dist)
  step_exposed_ei(sim$ppl, inf_dist = sim$inf_dist)
  step_transmission_se(sim$ppl, sim$nd, beta = beta, exp_dist = sim$exp_dist)
}

n_in_state <- function(sim, code) sum(sim$ppl$state == code)

# ── one SEIR run over a fixed horizon, returning the S / E / I / R trajectory ────
run_seir_traj <- function(n, n_seed, beta, exp_duration, inf_duration, nticks) {
  sim <- new_seir(n, n_seed, exp_duration, inf_duration)
  traj <- matrix(0L, nrow = nticks + 1L, ncol = 4L,
                 dimnames = list(NULL, c("S", "E", "I", "R")))
  traj[1L, ] <- c(n - n_seed, 0L, n_seed, 0L)
  for (tick in seq_len(nticks)) {
    seir_step(sim, beta)
    traj[tick + 1L, ] <- c(n_in_state(sim, states[["S"]]),
                           n_in_state(sim, states[["E"]]),
                           n_in_state(sim, states[["I"]]),
                           n_in_state(sim, states[["R"]]))
  }
  traj
}

# ── run to completion (no exposed or infectious agents left) -> attack fraction ──
final_attack_fraction <- function(n, n_seed, beta, exp_duration, inf_duration,
                                   max_ticks = 5000L) {
  sim <- new_seir(n, n_seed, exp_duration, inf_duration)
  for (tick in seq_len(max_ticks)) {
    seir_step(sim, beta)
    if (n_in_state(sim, states[["E"]]) + n_in_state(sim, states[["I"]]) == 0L) break
  }
  (n - n_in_state(sim, states[["S"]])) / n   # fraction ever infected
}

# ── Kermack–McKendrick final size ──────────────────────────────────────────────
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

N  <- 100000L
De <- 5L    # incubation (exposed) period in ticks
D  <- 10L   # infectious period in ticks

# ── 1. single-run trajectories ──────────────────────────────────────────────────
beta_traj <- 0.25
R0_traj   <- beta_traj * D            # effective R0 for SEIR (full period; see header)
nticks_traj <- 250L
traj_time <- system.time(
  traj <- run_seir_traj(N, n_seed = 100L, beta = beta_traj,
                        exp_duration = De, inf_duration = D, nticks = nticks_traj)
)
cat(sprintf("run_seir_traj: %s agents x %d ticks took %.3f s (%.1f ns/agent-tick)\n",
            format(N, big.mark = ","), nticks_traj, traj_time[["elapsed"]],
            1e9 * traj_time[["elapsed"]] / (as.numeric(N) * nticks_traj)))

png(file.path(out_dir, "seir_trajectories.png"), width = 1000, height = 650, res = 110)
ticks <- 0:(nrow(traj) - 1L)
matplot(ticks, traj, type = "l", lwd = 2.5, lty = 1,
        col = c("#2c7fb8", "#e6a000", "#d7301f", "#238b45"),
        xlab = "tick", ylab = "number of agents",
        main = sprintf("razer SEIR trajectories  (N = %s, R0 = beta*D = %.2f)",
                       format(N, big.mark = ","), R0_traj))
legend("right", legend = c("S", "E", "I", "R"), bty = "n", lwd = 2.5,
       col = c("#2c7fb8", "#e6a000", "#d7301f", "#238b45"))
dev.off()

# ── 2. attack-fraction sweep vs Kermack–McKendrick ───────────────────────────────
R0_grid <- seq(0.5, 4.0, by = 0.25)
sweep_time <- system.time(
  sim_af <- vapply(R0_grid, function(R0)
    final_attack_fraction(N, n_seed = 50L, beta = R0 / D,
                          exp_duration = De, inf_duration = D),
    numeric(1))
)
cat(sprintf("attack-fraction sweep: %d runs to completion (N = %s each) took %.3f s\n",
            length(R0_grid), format(N, big.mark = ","), sweep_time[["elapsed"]]))
km_af <- vapply(R0_grid, km_attack_fraction, numeric(1))

png(file.path(out_dir, "seir_attack_fraction.png"), width = 1000, height = 650, res = 110)
plot(R0_grid, km_af, type = "l", lwd = 2.5, col = "black",
     ylim = c(0, 1), xlab = expression(R[0]), ylab = "attack fraction",
     main = "Attack fraction: razer SEIR vs Kermack–McKendrick")
points(R0_grid, sim_af, pch = 19, cex = 1.1, col = "#d7301f")
abline(v = 1, lty = 3, col = "grey50")
legend("bottomright", bty = "n",
       legend = c("Kermack-McKendrick  A = 1 - exp(-R0 A)", "razer SEIR simulation"),
       lwd = c(2.5, NA), pch = c(NA, 19), col = c("black", "#d7301f"))
dev.off()

# ── comparison table ─────────────────────────────────────────────────────────────
cmp <- data.frame(
  R0                 = R0_grid,
  simulated          = round(sim_af, 4),
  kermack_mckendrick = round(km_af, 4),
  abs_error          = round(abs(sim_af - km_af), 4)
)
cat("\nAttack fraction: simulated SEIR vs Kermack-McKendrick\n")
print(cmp, row.names = FALSE)
cat(sprintf("\nMax absolute error for R0 >= 1.5: %.4f\n",
            max(cmp$abs_error[cmp$R0 >= 1.5])))
cat(sprintf("Plots written to: %s\n", normalizePath(out_dir)))
