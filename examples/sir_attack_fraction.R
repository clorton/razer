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

# `library(pkg)` attaches a package so its exported names resolve unqualified
# (like a wildcard `import`); the package name is given bare, not as a string.
library(razer)

# Seed R's global RNG so the stochastic run is reproducible. `L` makes an integer
# literal (a plain `20260604` would be a double).
set.seed(20260604L)

# `<-` is assignment. `laser_states()` returns a *named* integer vector; the
# comment shows its contents. Name-index it below with `states[["S"]]`.
states <- laser_states()   # c(S = 0, E = 1, I = 2, R = 3, D = -1)

# `function(args) body` is a first-class closure value, bound here to `new_sir`.
# ── build + seed a single-node SIR population ───────────────────────────────────
new_sir <- function(n, n_seed, inf_duration) {
  # `$new(...)` is the R6/extendr constructor call; `$` reaches the method.
  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", states[["S"]])
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)

  # Seed the first n_seed agents as infectious with a full infectious timer.
  # `rep(x, n)` makes a length-n vector of copies of x (like a fill).
  state0 <- rep(states[["S"]], n)
  # `seq_len(k)` is the integer vector 1..k (1-based); using it to index `state0`
  # selects the first k elements, and `<-` writes the scalar into all of them.
  state0[seq_len(n_seed)] <- states[["I"]]
  ppl$state <- state0   # `$<-` writes the whole column back into the frame
  timer0 <- rep(0L, n)
  timer0[seq_len(n_seed)] <- inf_duration
  ppl$timer <- timer0

  # Single well-mixed node holding the total population and an I tally.
  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  # `list(name = value, ...)` builds a heterogeneous, named record (like a struct
  # / dict); it is the function's return value (last expression returns implicitly).
  list(ppl = ppl, nd = nd,
       inf_dist = dist_constant(inf_duration),  # fixed infectious period
       imm_dist = dist_constant(0))             # SIR: no waning, R timer unused
}

# Advance one tick. Transitions run downstream-first: I->R recovery before S->I
# transmission, so an agent infected this tick is not decremented by the recovery
# step in the same tick (see header note on effective infectious period).
sir_step <- function(sim, beta) {
  # `sim$ppl` reads the `ppl` element of the list; `name = value` at a call site
  # passes a named argument (order-independent, like a keyword argument).
  step_infectious_ir(sim$ppl, imm_dist = sim$imm_dist)
  step_transmission_si(sim$ppl, sim$nd, beta = beta, inf_dist = sim$inf_dist)
}

# One-line closure. `sim$ppl$state == code` is a vectorized elementwise compare
# yielding a logical vector; `sum(...)` counts the TRUEs (TRUE coerces to 1).
n_in_state <- function(sim, code) sum(sim$ppl$state == code)

# ── one SIR run over a fixed horizon, returning the S / I / R trajectory ─────────
run_sir_traj <- function(n, n_seed, beta, inf_duration, nticks) {
  sim <- new_sir(n, n_seed, inf_duration)
  # `matrix(fill, nrow=, ncol=, dimnames=)` allocates a 2-D array (stored
  # COLUMN-MAJOR internally). `dimnames = list(rownames, colnames)`; `NULL` leaves
  # rows unnamed, and the column-name vector lets us see S/I/R headers.
  traj <- matrix(0L, nrow = nticks + 1L, ncol = 3L,
                 dimnames = list(NULL, c("S", "I", "R")))
  # `traj[i, ]` selects row i, all columns; assigning a length-3 vector fills it.
  traj[1L, ] <- c(n - n_seed, n_seed, 0L)
  for (tick in seq_len(nticks)) {           # `for (x in vec)` foreach loop
    sir_step(sim, beta)
    traj[tick + 1L, ] <- c(n_in_state(sim, states[["S"]]),
                           n_in_state(sim, states[["I"]]),
                           n_in_state(sim, states[["R"]]))
  }
  traj                                       # implicit return of the matrix
}

# ── run to completion (no infectious agents left) and return the attack fraction ─
final_attack_fraction <- function(n, n_seed, beta, inf_duration, max_ticks = 5000L) {
  sim <- new_sir(n, n_seed, inf_duration)
  for (tick in seq_len(max_ticks)) {
    sir_step(sim, beta)
    if (n_in_state(sim, states[["I"]]) == 0L) break   # `break` exits the loop
  }
  # Parenthesized expression is the implicit return value.
  (n - n_in_state(sim, states[["S"]])) / n   # fraction ever infected
}

# ── Kermack–McKendrick final size ──────────────────────────────────────────────
# Solve A = 1 - exp(-R0 * A) for the non-trivial root in (0, 1). For R0 <= 1 the
# only root is A = 0 (no epidemic).
km_attack_fraction <- function(R0) {
  if (R0 <= 1) return(0)
  # Nested closure capturing `R0` (lexical scope); `f(A)` is the residual to zero.
  f <- function(A) 1 - exp(-R0 * A) - A
  # `uniroot` is base R's 1-D root finder on the bracket [lower, upper]; it returns
  # a list, and `$root` pulls out the located root.
  uniroot(f, lower = 1e-9, upper = 1 - 1e-12)$root
}

# ── output location (next to this script) ───────────────────────────────────────
# Idiom for "find the directory of the running script". `commandArgs(FALSE)`
# returns ALL launch args including Rscript's own `--file=<path>`.
args       <- commandArgs(trailingOnly = FALSE)
# `grep(pat, x, value = TRUE)` returns the matching *elements* (not indices).
file_arg   <- grep("^--file=", args, value = TRUE)
# `if (...) a else b` is an expression that yields a value. `length(x)` is truthy
# when non-zero. `sub(pat, repl, x)` strips the `--file=` prefix; `dirname` takes
# the directory part. The `else "."` fallback covers being sourced interactively.
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")   # OS-correct path join
# `dir.create` mkdir; `recursive` = mkdir -p, `showWarnings = FALSE` = no warn if it exists.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N <- 100000L
D <- 10L

# ── 1. single-run trajectories ──────────────────────────────────────────────────
beta_traj <- 0.25
R0_traj   <- beta_traj * (D - 1)        # effective R0 (see header)
nticks_traj <- 200L
# `system.time(expr)` runs `expr` and returns a named numeric vector of timings;
# the `traj <- ...` assignment inside still takes effect in this scope.
traj_time <- system.time(
  traj <- run_sir_traj(N, n_seed = 100L, beta = beta_traj, inf_duration = D, nticks = nticks_traj)
)
# `cat` writes to stdout (no quoting/newline added). `sprintf` is C-style format:
# %s string, %d integer, %.3f fixed-point. `format(N, big.mark=",")` adds thousands
# separators. `traj_time[["elapsed"]]` name-indexes the wall-clock field of the vector.
cat(sprintf("run_sir_traj: %s agents x %d ticks took %.3f s (%.1f ns/agent-tick)\n",
            format(N, big.mark = ","), nticks_traj, traj_time[["elapsed"]],
            # `as.numeric` casts the integer N to double to avoid integer overflow.
            1e9 * traj_time[["elapsed"]] / (as.numeric(N) * nticks_traj)))

# Graphics-device pattern: open a device (here a PNG file), draw into it with the
# plot calls below, then close it with `dev.off()` to flush the file.
png(file.path(out_dir, "sir_trajectories.png"), width = 1000, height = 650, res = 110)
# `a:b` is the inclusive integer-range operator; `nrow` is the matrix's row count.
ticks <- 0:(nrow(traj) - 1L)
# `matplot` plots each COLUMN of the matrix as its own line vs the x vector.
matplot(ticks, traj, type = "l", lwd = 2.5, lty = 1,
        col = c("#2c7fb8", "#d7301f", "#238b45"),
        xlab = "tick", ylab = "number of agents",
        main = sprintf("razer SIR trajectories  (N = %s, R0 = beta*(D-1) = %.2f)",
                       format(N, big.mark = ","), R0_traj))
# `legend` overlays a key; `bty = "n"` draws no box around it.
legend("right", legend = c("S", "I", "R"), bty = "n", lwd = 2.5,
       col = c("#2c7fb8", "#d7301f", "#238b45"))
dev.off()   # close the device -> writes the PNG

# ── 2. attack-fraction sweep vs Kermack–McKendrick ───────────────────────────────
# `seq(from, to, by =)` builds an arithmetic sequence (the R0 grid).
R0_grid <- seq(0.5, 4.0, by = 0.25)
sweep_time <- system.time(
  # `vapply(x, f, template)` is a typed map: apply f to each element, asserting each
  # result matches `numeric(1)` (one double) -> returns a numeric vector. The inline
  # `function(R0) ...` is the per-element callback (its `{}` are optional here).
  sim_af <- vapply(R0_grid, function(R0)
    final_attack_fraction(N, n_seed = 50L, beta = R0 / (D - 1), inf_duration = D),
    numeric(1))
)
cat(sprintf("attack-fraction sweep: %d runs to completion (N = %s each) took %.3f s\n",
            length(R0_grid), format(N, big.mark = ","), sweep_time[["elapsed"]]))
# Passing the named function directly (no wrapper lambda needed when arity matches).
km_af   <- vapply(R0_grid, km_attack_fraction, numeric(1))

png(file.path(out_dir, "sir_attack_fraction.png"), width = 1000, height = 650, res = 110)
# `expression(R[0])` is a plotmath label: renders as R with subscript 0.
plot(R0_grid, km_af, type = "l", lwd = 2.5, col = "black",
     ylim = c(0, 1), xlab = expression(R[0]), ylab = "attack fraction",
     main = "Attack fraction: razer SIR vs Kermack–McKendrick")
# `points` overlays markers on the existing plot; `abline(v = 1)` adds a vertical line.
points(R0_grid, sim_af, pch = 19, cex = 1.1, col = "#d7301f")
abline(v = 1, lty = 3, col = "grey50")
legend("bottomright", bty = "n",
       legend = c("Kermack-McKendrick  A = 1 - exp(-R0 A)", "razer simulation"),
       lwd = c(2.5, NA), pch = c(NA, 19), col = c("black", "#d7301f"))
dev.off()

# ── comparison table ─────────────────────────────────────────────────────────────
# `data.frame` is R's tabular type (named columns of equal length, like a SQL row
# set / pandas DataFrame). `round(x, 4)` and `abs(...)` are vectorized over columns.
cmp <- data.frame(
  R0                 = R0_grid,
  simulated          = round(sim_af, 4),
  kermack_mckendrick = round(km_af, 4),
  abs_error          = round(abs(sim_af - km_af), 4)
)
cat("\nAttack fraction: simulated vs Kermack-McKendrick\n")
print(cmp, row.names = FALSE)   # suppress the auto-generated 1..n row labels
# `cmp$col` selects a column by name. `cmp$R0 >= 1.5` is a logical vector;
# `cmp$abs_error[<logical>]` keeps the elements where it is TRUE (boolean masking).
cat(sprintf("\nMax absolute error for R0 >= 1.5: %.4f\n",
            max(cmp$abs_error[cmp$R0 >= 1.5])))
# `normalizePath` resolves to an absolute, canonical path for the message.
cat(sprintf("Plots written to: %s\n", normalizePath(out_dir)))
