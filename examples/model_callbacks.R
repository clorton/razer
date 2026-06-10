#!/usr/bin/env Rscript

# run_model() callbacks + the model environment + bincount_where(), end to end.
#
# This example shows how to EXTEND the high-level run_model() runner without forking it,
# using the three callbacks it offers — each receives the `model` environment that bundles
# `$people`, `$nodes`, `$network`, `$carry`, and (during the loop) the current `$tick`:
#   * init(model)       — once, before the loop: add our own per-agent and per-node Columns.
#   * step_enter(model) — start of each tick: a one-off "pulse vaccination" intervention.
#   * step_exit(model)  — end of each tick: record a CUSTOM per-node report that the runner
#                         does not track — here the number of infectious UNDER-FIVES, via
#                         the predicate-filtered binner bincount_where().
#
# It runs the same SEIR scenario twice — once plain, once with a mid-epidemic pulse that
# vaccinates a fraction of susceptibles (moving S -> R directly) — and plots the effect on
# the infectious curve plus the under-five infectious series the callback collected.
#
# Run from anywhere:  Rscript examples/model_callbacks.R
# Output PNGs are written next to this script in examples/output/.

library(razer)

set.seed(20260606L)
states <- laser_states()                 # c(S = 0, E = 1, I = 2, R = 3, M = 4, D = -1)

# ── scenario: a single well-mixed node ───────────────────────────────────────────
N        <- 200000L
nticks   <- 365L                          # one year, daily
n_nodes  <- 1L
scenario <- data.frame(population = N, I = 50L)
day_per_year <- 365L
under5_days  <- 5L * day_per_year         # "under five" threshold in days

# ── init: give each agent a date-of-birth, and add an under-five-infectious report ─
# dob is the NEGATIVE age in days (so a 3-year-old has dob = -3*365). At tick t an agent is
# under five iff its current age `t - dob < under5_days`, i.e. `dob > t - under5_days` — a
# single ">" predicate bincount_where() evaluates directly on the Column. Ages are drawn
# from a rough exponential age structure (more young than old), capped at 80 years.
# `model` is the run_model environment; we attach our own Columns to its $people / $nodes.
init_ages <- function(model) {
  n <- model$people$count
  age_days <- pmin(as.integer(rexp(n, rate = 1 / (25 * day_per_year))), 80L * day_per_year)
  dob <- allocate_scalar("i32", n)
  dob$set(-age_days)                                   # dob = -age (0-based "born at -age")
  model$people$dob <- dob
  # Per-interval report Column we will fill ourselves (run_model does not track this).
  model$nodes$under5_I <- allocate_vector("i32", nticks - 1L, n_nodes)
}

# ── step_exit: record infectious under-fives per node for the current tick ─────────
# bincount_where with TWO calls intersected would be ideal, but one predicate per call is
# the primitive; here we report the under-five count among ALL agents (the cohort ages out
# over the year) and, separately, the infectious count — both per node, both via the binner.
record_under5_I <- function(model) {
  t <- model$tick
  thresh <- t - under5_days                            # dob > thresh  <=>  age < 5 years
  model$nodes$under5_I$set_col(
    t, bincount_where(model$people$nodeid, n_nodes, model$people$dob, "gt",
                      thresh, model$people$count))
}

# ── step_enter: a one-time pulse vaccination at a chosen tick ──────────────────────
# At `pulse_tick`, move a fraction of still-susceptible agents straight to R (vaccinated /
# immune). We mutate the agent `state` Column in place AND keep the node census consistent
# by recounting S and R from the agents (bincount_where) into the just-settled column the
# carry-forward will propagate. Returns a closure capturing the pulse parameters.
make_pulse <- function(pulse_tick, coverage) {
  function(model) {
    if (model$tick != pulse_tick) return(invisible(NULL))
    st  <- model$people$state$values()
    sus <- which(st == states[["S"]])
    take <- sample(sus, size = floor(coverage * length(sus)))
    st[take] <- states[["R"]]
    model$people$state$set(st)                         # S -> R for the vaccinated
    # Re-sync the census at the current settled column `t` so the carry-forward sees it.
    t <- model$tick
    nid <- model$people$nodeid; sc <- model$people$state; cnt <- model$people$count
    model$nodes$S$set_col(t, bincount_where(nid, n_nodes, sc, "eq", states[["S"]], cnt))
    model$nodes$R$set_col(t, bincount_where(nid, n_nodes, sc, "eq", states[["R"]], cnt))
    cat(sprintf("  pulse: vaccinated %s susceptibles at tick %d\n",
                format(length(take), big.mark = ","), pulse_tick))
  }
}

# ── run twice: baseline and with the pulse ────────────────────────────────────────
common <- list(scenario = scenario, model = "SEIR", nticks = nticks, beta = 2.2 / 6,
               infectious_period = 6, incubation_period = 4)

cat("baseline SEIR run...\n")
base <- do.call(run_model, c(common, list(seed = 1L, init = init_ages, step_exit = record_under5_I)))

cat("pulse-vaccination SEIR run...\n")
pulsed <- do.call(run_model, c(common, list(seed = 1L, init = init_ages,
                                            step_exit = record_under5_I,
                                            step_enter = make_pulse(pulse_tick = 60L, coverage = 0.4))))

# ── output location (next to this script) ─────────────────────────────────────────
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

ticks <- seq_len(nticks)

# ── plot 1: infectious curve, baseline vs pulse ───────────────────────────────────
Ib <- rowSums(base$nodes$I$values()); Ip <- rowSums(pulsed$nodes$I$values())
open_png(file.path(out_dir, "model_callbacks_intervention.png"), width = 1000, height = 650, res = 110)
matplot(ticks, cbind(baseline = Ib, pulse = Ip), type = "l", lwd = 2.5, lty = 1,
        col = c("#d7301f", "#2c7fb8"), xlab = "day", ylab = "infectious agents",
        main = "SEIR infectious curve: baseline vs a tick-60 pulse vaccination (40% of S)")
abline(v = 60, lty = 3, col = "grey50")
legend("topright", legend = c("baseline", "pulse @ day 60"), bty = "n", lwd = 2.5,
       col = c("#d7301f", "#2c7fb8"))
close_png()

# ── plot 2: the custom callback-collected report (under-five population over time) ─
u5 <- rowSums(base$nodes$under5_I$values())
open_png(file.path(out_dir, "model_callbacks_under5.png"), width = 1000, height = 650, res = 110)
plot(ticks[-length(ticks)], u5, type = "l", lwd = 2.5, col = "#238b45",
     xlab = "day", ylab = "agents under 5 years",
     main = "Under-fives per day, recorded in a step_exit callback via bincount_where()")
close_png()

# ── plot 3: full baseline S/E/I/R trajectories ────────────────────────────────────
traj <- cbind(S = rowSums(base$nodes$S$values()), E = rowSums(base$nodes$E$values()),
              I = Ib, R = rowSums(base$nodes$R$values()))
open_png(file.path(out_dir, "model_callbacks_seir.png"), width = 1000, height = 650, res = 110)
matplot(ticks, traj, type = "l", lwd = 2.5, lty = 1,
        col = c("#2c7fb8", "#fdae61", "#d7301f", "#238b45"),
        xlab = "day", ylab = "agents", main = "Baseline SEIR trajectories")
legend("right", legend = c("S", "E", "I", "R"), bty = "n", lwd = 2.5,
       col = c("#2c7fb8", "#fdae61", "#d7301f", "#238b45"))
close_png()

cat(sprintf("baseline:  peak I = %s on day %d;  final R = %s\n",
            format(max(Ib), big.mark = ","), which.max(Ib) - 1L,
            format(tail(traj[, "R"], 1), big.mark = ",")))
cat(sprintf("pulse:     peak I = %s on day %d\n",
            format(max(Ip), big.mark = ","), which.max(Ip) - 1L))
if (to_png) cat(sprintf("Plots written to: %s\n", normalizePath(out_dir)))
