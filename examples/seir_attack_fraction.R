#!/usr/bin/env Rscript

# SEIR example: trajectories and the Kermack–McKendrick final-size relation.
#
# Companion to sir_attack_fraction.R. It
#   1. runs razer's agent-based SEIR dynamics for a single parameter set and plots
#      the S / E / I / R trajectories over time, then
#   2. sweeps the basic reproduction number R0 and compares the simulated attack
#      fraction (fraction of the population ever infected) against the
#      Kermack–McKendrick final-size equation
#
#          A = 1 - exp(-R0 * A).
#
# The latent (exposed) period only *delays* the epidemic; it does not change the
# final size, so the SEIR attack fraction obeys the same relation as SIR with R0
# determined by the infectious period alone.
#
# Built on the Column kernels. SEIR is `step_sir` with the absorbing state set to R:
# it performs E->I (drawing the infectious timer) and I->R, returning per-node counts
# the model applies to its E/I/R census; `transmission` performs S->E (drawing the
# incubation timer) and returns new infections per node. Per-tick order is the single
# razer ordering: carry_forward -> step_sir -> calc_foi -> transmission, with `calc_foi`
# immediately before `transmission`; it reads the settled start-of-interval census I[t],
# so each infectious agent contributes on exactly its D census columns and
#
#     R0 = beta * D.
#
# Run from anywhere:  Rscript examples/seir_attack_fraction.R
# Output PNGs are written next to this script in examples/output/.

library(razer)

set.seed(20260604L)

states <- laser_states()   # c(S = 0, E = 1, I = 2, R = 3, M = 4, D = -1)

# ── one single-node SEIR run over `horizon` recorded states ──────────────────────────
# Column-based well-mixed population (u16 timer). SEIR is `step_sir` with absorbing = R;
# there is no M compartment, so its `waned` count is ignored. `stop_extinct = TRUE`
# breaks once E + I == 0.
run_seir_columns <- function(n, n_seed, exp_duration, inf_duration, beta_val, horizon,
                             stop_extinct = FALSE) {
  nn <- 1L
  state  <- allocate_scalar("u8",  n)
  timer  <- allocate_scalar("u16", n)        # incubation/infectious timer (u16)
  nodeid <- allocate_scalar("u16", n)
  s0 <- rep(states[["S"]], n); s0[seq_len(n_seed)] <- states[["I"]]; state$set(s0)
  tm <- rep(0L, n);            tm[seq_len(n_seed)] <- inf_duration; timer$set(tm)

  # No M (maternal) compartment: step_sir's `waned` (M->S) count is always 0 here, so we
  # simply never allocate or update an M census — the return-counts kernels let a model
  # carry only the compartments it actually has.
  S <- allocate_vector("i32", horizon, nn)
  E <- allocate_vector("i32", horizon, nn)
  I <- allocate_vector("i32", horizon, nn)
  R <- allocate_vector("i32", horizon, nn)
  foi <- allocate_vector("f64", horizon - 1L, nn)
  zeros <- rep(0L, (horizon - 1L) * nn)
  S$set(c(n - n_seed, zeros)); E$set(c(0L, zeros)); I$set(c(n_seed, zeros)); R$set(c(0L, zeros))
  N      <- values_map(n,        horizon, nn)
  beta   <- values_map(beta_val, horizon, nn)
  season <- values_map(1,        horizon, nn)
  net    <- matrix(0, nn, nn)
  exp_dist <- dist_constant(exp_duration)    # incubation period (S->E timer)
  inf_dist <- dist_constant(inf_duration)    # infectious period (E->I timer)

  last <- horizon - 1L
  for (tick in seq_len(horizon - 1L)) {
    t <- tick - 1L
    carry_forward_states(list(S, E, I, R), t)
    # step_sir(absorbing = R): E->I onset and I->R recovery. (waned is 0 — no M.)
    prog <- step_sir(state, timer, nodeid, n, nn, inf_dist, states[["R"]])
    move_count(E, I, prog$onset,   t)
    move_count(I, R, prog$cleared, t)
    calc_foi(I, N, beta, season, net, foi, t)                       # reads settled I[t]; just before transmit
    inf <- transmission(state, timer, nodeid, n, foi, t, states[["E"]], exp_dist)
    move_count(S, E, inf, t)
    if (stop_extinct && sum(E$col(tick)) + sum(I$col(tick)) == 0L) { last <- tick; break }
  }
  list(S = S, E = E, I = I, R = R, last = last)
}

# ── one SEIR run over a fixed horizon, returning the S / E / I / R trajectory matrix ─
run_seir_traj <- function(n, n_seed, beta, exp_duration, inf_duration, nticks) {
  res <- run_seir_columns(n, n_seed, exp_duration, inf_duration, beta, nticks + 1L)
  cbind(S = res$S$values()[, 1], E = res$E$values()[, 1],
        I = res$I$values()[, 1], R = res$R$values()[, 1])
}

# ── run to completion (no exposed or infectious agents left) -> attack fraction ──────
final_attack_fraction <- function(n, n_seed, beta, exp_duration, inf_duration,
                                   max_ticks = 5000L) {
  res <- run_seir_columns(n, n_seed, exp_duration, inf_duration, beta, max_ticks + 1L,
                          stop_extinct = TRUE)
  (n - sum(res$S$col(res$last))) / n
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
beta_traj   <- 0.25
R0_traj     <- beta_traj * D            # effective R0 for SEIR (full period; see header)
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
