#!/usr/bin/env Rscript

# SIR example: trajectories and the Kermack–McKendrick final-size relation.
#
# This script
#   1. runs razer's agent-based SIR kernels for a single parameter set and plots
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
# Built on the Column kernels (`calc_foi` / `step_sir` / `transmission`). Per-tick order is
# the single razer ordering: carry_forward -> step_sir (I->R) -> calc_foi -> transmission
# (S->I), with `calc_foi` placed IMMEDIATELY before `transmission` (no step between them).
# `calc_foi` reads the SETTLED start-of-interval infectious census I[t], so each infectious
# agent contributes to the force of infection on exactly the D census columns it occupies
# (it enters I one column after infection and recovers D columns later). The effective
# basic reproduction number is therefore the full infectious period,
#
#     R0 = beta * D
#
# (with the per-tick force of infection p = 1 - exp(-beta*I/N)). With this mapping the
# deterministic final size satisfies the equation above, run to completion — which the
# sweep below does by iterating until no infectious agents remain.
#
# Run from anywhere:  Rscript examples/sir_attack_fraction.R
# Output PNGs are written next to this script in examples/output/.

# `library(pkg)` attaches a package so its exported names resolve unqualified
# (like a wildcard `import`); the package name is given bare, not as a string.
library(razer)

# Seed R's global RNG (the within-tick draws use a separate, non-seedable RNG, but
# this keeps the seeding/jitter reproducible). `L` makes an integer literal.
set.seed(20260604L)

# `<-` is assignment. `laser_states()` returns a *named* integer vector; name-index it
# below with `states[["S"]]`.
states <- laser_states()   # c(S = 0, E = 1, I = 2, R = 3, M = 4, D = -1)

# ── one single-node SIR run over `horizon` recorded states ───────────────────────────
# Builds a Column-based well-mixed population and advances it. The S/I/R census is
# maintained incrementally in `horizon`-column buffers. `stop_extinct = TRUE` breaks as
# soon as the infectious count hits zero (for the attack-fraction sweep); the trajectory
# run leaves it FALSE so post-epidemic columns hold the steady (flat) S/R values.
# Returns the census Columns and the last completed column index.
run_sir_columns <- function(n, n_seed, beta_val, D, horizon, stop_extinct = FALSE) {
  nn <- 1L                                   # single well-mixed node
  # Per-agent arrays (state u8, timer u16, node id u16 — all node 0).
  state  <- allocate_scalar("u8",  n)
  timer  <- allocate_scalar("u16", n)
  nodeid <- allocate_scalar("u16", n)
  s0 <- rep(states[["S"]], n); s0[seq_len(n_seed)] <- states[["I"]]; state$set(s0)
  tm <- rep(0L, n);            tm[seq_len(n_seed)] <- D;            timer$set(tm)

  # Census (horizon x 1) only — the kernels return counts; no flow buffers needed here.
  S <- allocate_vector("i32", horizon, nn)
  I <- allocate_vector("i32", horizon, nn)
  R <- allocate_vector("i32", horizon, nn)
  foi <- allocate_vector("f64", horizon - 1L, nn)
  zeros <- rep(0L, (horizon - 1L) * nn)
  S$set(c(n - n_seed, zeros)); I$set(c(n_seed, zeros)); R$set(c(0L, zeros))
  # Constant population (no vital dynamics) and a flat transmission grid; values_map
  # builds the n_ticks x n_nodes grids calc_foi reads.
  N      <- values_map(n,        horizon, nn)
  beta   <- values_map(beta_val, horizon, nn)
  season <- values_map(1,        horizon, nn)
  net    <- matrix(0, nn, nn)                # single node: no spatial coupling
  inf_dist <- dist_constant(D)               # fixed infectious period

  last <- horizon - 1L
  for (tick in seq_len(horizon - 1L)) {
    t <- tick - 1L
    carry_forward_states(list(S, I, R), t)              # copy column t -> t+1
    # step_sir with absorbing = R: I->R recovery (no M/E agents, so waned/onset are 0).
    rec <- step_sir(state, timer, nodeid, n, nn, inf_dist, states[["R"]])
    move_count(I, R, rec$cleared, t)
    calc_foi(I, N, beta, season, net, foi, t)           # reads settled I[t]; just before transmit -> R0 = beta*D
    # transmission S->I returns new infections per node; apply the S->I census delta.
    inf <- transmission(state, timer, nodeid, n, foi, t, states[["I"]], inf_dist)
    move_count(S, I, inf, t)
    if (stop_extinct && sum(I$col(tick)) == 0L) { last <- tick; break }  # I[tick] is current
  }
  list(S = S, I = I, R = R, last = last)
}

# ── one SIR run over a fixed horizon, returning the S / I / R trajectory matrix ──────
run_sir_traj <- function(n, n_seed, beta, inf_duration, nticks) {
  res <- run_sir_columns(n, n_seed, beta, inf_duration, nticks + 1L)
  # `$values()` copies each census buffer back as an (nticks+1) x 1 matrix; `[, 1]` is
  # its single node column. `cbind` glues the three into one matrix with named columns.
  cbind(S = res$S$values()[, 1], I = res$I$values()[, 1], R = res$R$values()[, 1])
}

# ── run to completion (no infectious agents left) and return the attack fraction ─────
final_attack_fraction <- function(n, n_seed, beta, inf_duration, max_ticks = 5000L) {
  res <- run_sir_columns(n, n_seed, beta, inf_duration, max_ticks + 1L, stop_extinct = TRUE)
  # Fraction ever infected = 1 - (final susceptibles / N). `$col(last)` reads the final
  # recorded S column (a length-1 node vector here).
  (n - sum(res$S$col(res$last))) / n
}

# ── Kermack–McKendrick final size ──────────────────────────────────────────────
# Solve A = 1 - exp(-R0 * A) for the non-trivial root in (0, 1). For R0 <= 1 the
# only root is A = 0 (no epidemic).
km_attack_fraction <- function(R0) {
  if (R0 <= 1) return(0)
  f <- function(A) 1 - exp(-R0 * A) - A             # residual to zero (captures R0)
  uniroot(f, lower = 1e-9, upper = 1 - 1e-12)$root  # base R 1-D root finder
}

# ── output location (next to this script) ───────────────────────────────────────
args       <- commandArgs(trailingOnly = FALSE)        # all launch args incl. --file=
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N <- 100000L
D <- 10L

# ── 1. single-run trajectories ──────────────────────────────────────────────────
beta_traj   <- 0.25
R0_traj     <- beta_traj * D            # effective R0 (full infectious period; see header)
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
        main = sprintf("razer SIR trajectories  (N = %s, R0 = beta*D = %.2f)",
                       format(N, big.mark = ","), R0_traj))
legend("right", legend = c("S", "I", "R"), bty = "n", lwd = 2.5,
       col = c("#2c7fb8", "#d7301f", "#238b45"))
dev.off()

# ── 2. attack-fraction sweep vs Kermack–McKendrick ───────────────────────────────
R0_grid <- seq(0.5, 4.0, by = 0.25)
sweep_time <- system.time(
  sim_af <- vapply(R0_grid, function(R0)
    final_attack_fraction(N, n_seed = 50L, beta = R0 / D, inf_duration = D),
    numeric(1))
)
cat(sprintf("attack-fraction sweep: %d runs to completion (N = %s each) took %.3f s\n",
            length(R0_grid), format(N, big.mark = ","), sweep_time[["elapsed"]]))
km_af <- vapply(R0_grid, km_attack_fraction, numeric(1))

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
