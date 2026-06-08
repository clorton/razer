#!/usr/bin/env Rscript

# SEIR example: trajectories and the Kermack–McKendrick final-size relation.
#
# Companion to sir_attack_fraction.R. It
#   1. runs razer's SEIR model with the high-level runner run_model() for a single
#      parameter set and plots the S / E / I / R trajectories over time, then
#   2. sweeps the basic reproduction number R0 and compares the simulated attack
#      fraction (fraction of the population ever infected) against the
#      Kermack–McKendrick final-size equation
#
#          A = 1 - exp(-R0 * A).
#
# The latent (exposed) period only *delays* the epidemic; it does not change the final
# size, so the SEIR attack fraction obeys the same relation as SIR with R0 determined by
# the infectious period alone. run_model() advances the Column kernels in the single razer
# per-tick order (carry_forward -> step -> calc_foi -> transmission); because calc_foi
# reads the settled start-of-interval infectious census, each infectious agent contributes
# to the force of infection on exactly its D census columns, so
#
#     R0 = beta * D
#
# (run_model derives beta = R0 / mean(infectious_period)). See simple_sir.R /
# endemic_sir.R / engwal_measles.R for hand-wired loops that go beyond run_model's
# closed-population menagerie (vital dynamics, importation, a maternal compartment).
#
# Run from anywhere:  Rscript examples/seir_attack_fraction.R
# Output PNGs are written next to this script in examples/output/.

library(razer)

N  <- 100000L
De <- 5L    # incubation (exposed) period in ticks
D  <- 10L   # infectious period in ticks

# ── one SEIR run via run_model(), returning the S / E / I / R trajectory matrix ──────
# SEIR adds the latent E stage, so run_model needs an `incubation_period` in addition to
# the `infectious_period`; both are passed as bare numbers (coerced to constants). `seed`
# makes the stochastic run reproducible.
run_seir_traj <- function(n, n_seed, r0, De, D, nticks) {
  m <- run_model(data.frame(population = n, I = n_seed), "SEIR",
                 nticks = nticks, r0 = r0, infectious_period = D,
                 incubation_period = De, seed = 1L)
  cbind(S = m$nodes$S$values()[, 1], E = m$nodes$E$values()[, 1],
        I = m$nodes$I$values()[, 1], R = m$nodes$R$values()[, 1])
}

# ── attack fraction = fraction ever infected = 1 - final_S / N ───────────────────────
# Fixed horizon long enough for the (E-delayed) epidemic to complete.
run_to_completion_ticks <- 2000L
final_attack_fraction <- function(n, n_seed, r0, De, D, nticks = run_to_completion_ticks) {
  m <- run_model(data.frame(population = n, I = n_seed), "SEIR",
                 nticks = nticks, r0 = r0, infectious_period = D,
                 incubation_period = De, seed = 1L)
  S <- m$nodes$S$values()[, 1]
  (n - S[length(S)]) / n
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
# Device-aware: write a PNG when run non-interactively (Rscript); draw to the active
# device (e.g. RStudio's Plots pane) when sourced interactively.
to_png    <- !interactive()
open_png  <- function(path, ...) if (to_png) grDevices::png(path, ...)
close_png <- function() if (to_png) grDevices::dev.off()
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. single-run trajectories ──────────────────────────────────────────────────
R0_traj     <- 2.5                      # beta = R0 / D = 0.25
nticks_traj <- 250L
traj_time <- system.time(
  traj <- run_seir_traj(N, n_seed = 100L, r0 = R0_traj, De = De, D = D, nticks = nticks_traj)
)
cat(sprintf("run_seir_traj: %s agents x %d ticks took %.3f s (%.1f ns/agent-tick)\n",
            format(N, big.mark = ","), nticks_traj, traj_time[["elapsed"]],
            1e9 * traj_time[["elapsed"]] / (as.numeric(N) * nticks_traj)))

open_png(file.path(out_dir, "seir_trajectories.png"), width = 1000, height = 650, res = 110)
ticks <- 0:(nrow(traj) - 1L)
matplot(ticks, traj, type = "l", lwd = 2.5, lty = 1,
        col = c("#2c7fb8", "#e6a000", "#d7301f", "#238b45"),
        xlab = "tick", ylab = "number of agents",
        main = sprintf("razer SEIR trajectories  (N = %s, R0 = beta*D = %.2f)",
                       format(N, big.mark = ","), R0_traj))
legend("right", legend = c("S", "E", "I", "R"), bty = "n", lwd = 2.5,
       col = c("#2c7fb8", "#e6a000", "#d7301f", "#238b45"))
close_png()

# ── 2. attack-fraction sweep vs Kermack–McKendrick ───────────────────────────────
R0_grid <- seq(0.5, 4.0, by = 0.25)
sweep_time <- system.time(
  sim_af <- vapply(R0_grid, function(R0)
    final_attack_fraction(N, n_seed = 50L, r0 = R0, De = De, D = D),
    numeric(1))
)
cat(sprintf("attack-fraction sweep: %d run_model runs (N = %s, %d ticks each) took %.3f s\n",
            length(R0_grid), format(N, big.mark = ","), run_to_completion_ticks,
            sweep_time[["elapsed"]]))
km_af <- vapply(R0_grid, km_attack_fraction, numeric(1))

open_png(file.path(out_dir, "seir_attack_fraction.png"), width = 1000, height = 650, res = 110)
plot(R0_grid, km_af, type = "l", lwd = 2.5, col = "black",
     ylim = c(0, 1), xlab = expression(R[0]), ylab = "attack fraction",
     main = "Attack fraction: razer SEIR vs Kermack–McKendrick")
points(R0_grid, sim_af, pch = 19, cex = 1.1, col = "#d7301f")
abline(v = 1, lty = 3, col = "grey50")
legend("bottomright", bty = "n",
       legend = c("Kermack-McKendrick  A = 1 - exp(-R0 A)", "razer run_model() simulation"),
       lwd = c(2.5, NA), pch = c(NA, 19), col = c("black", "#d7301f"))
close_png()

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
if (to_png) cat(sprintf("Plots written to: %s\n", normalizePath(out_dir)))
