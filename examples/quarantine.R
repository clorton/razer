#!/usr/bin/env Rscript

# Quarantine / case isolation, via a user-defined "Q" agent state. A testing programme runs
# PERIODICALLY (every `test_period` ticks — here a fortnight) and is LEAKY: on a test day it
# detects each currently-infectious agent with probability `sensitivity` (so 1 - sensitivity
# are false negatives and missed), moving the detected ones I -> Q. We model no false
# positives. A quarantined agent is isolated — Q is not the I state, and the force of
# infection reads only the I census, so quarantined cases no longer transmit. After a fixed
# `quarantine_days` isolation period a step_update callback converts Q -> R via the generic
# step_timer_expire kernel.
#
# Q is registered with run_model(extra_states = "Q"); the disease kernels leave Q agents
# untouched, so both the isolation (I->Q) and the release (Q->R) are composed from callbacks
# with no kernel change. We run the same epidemic twice — without and with quarantine.
#
# Run from anywhere:  Rscript examples/quarantine.R
# In RStudio: source() this file and the plots appear in the Plots pane (no PNG is written);
# run via Rscript and it writes the PNG to examples/output/ instead (see `to_png` below).

library(razer)

set.seed(20260608L)

# ── parameters ────────────────────────────────────────────────────────────────────
N               <- 500000L
nticks          <- 200L
r0              <- 2.5
inf_duration    <- dist_gamma(2, 4)     # infectious period (mean 8 ticks)
test_period     <- 14L                   # test the population every fortnight (configurable)
sensitivity     <- 0.80                  # leaky detection: 80% of infectious caught per test
quarantine_days <- 10L                   # isolation period, then Q -> R
scenario        <- data.frame(population = N, I = 100L)
test_ticks      <- seq(0L, nticks - 2L, by = test_period)   # the fortnightly testing schedule

# ── intervention callbacks ──────────────────────────────────────────────────────────
# A small environment records the per-test-day detection yield for a plot (closures can
# write into an environment by reference).
track <- new.env(); track$day <- integer(0); track$detected <- integer(0)

# detect_and_isolate (step_exit): ONLY on test days, flip each infectious agent to Q with
# probability `sensitivity` (leaky), set its isolation timer, and apply the I -> Q delta.
detect_and_isolate <- function(model) {
  if (!(model$tick %in% test_ticks)) return(invisible(NULL))
  I <- model$states[["I"]]; Q <- model$states[["Q"]]
  cnt   <- model$people$count
  state <- model$people$state$values()
  node  <- model$people$nodeid$values()
  infectious <- which(state[seq_len(cnt)] == I)
  detected   <- infectious[stats::runif(length(infectious)) < sensitivity]   # imperfect test
  track$day      <- c(track$day, model$tick)
  track$detected <- c(track$detected, length(detected))
  if (!length(detected)) return(invisible(NULL))
  state[detected] <- Q
  timer <- model$people$timer$values()
  timer[detected] <- quarantine_days                       # isolation countdown -> Q->R
  model$people$state$set(state)
  model$people$timer$set(timer)
  counts <- tabulate(node[detected] + 1L, model$nodes$count)
  move_count(model$nodes$I, model$nodes$Q, counts, model$tick)   # I -> Q at slice tick+1
}

# release_recovered (step_update): Q -> R when the isolation timer expires (generic kernel).
release_recovered <- function(model) {
  Q <- model$states[["Q"]]; R <- model$states[["R"]]
  released <- step_timer_expire(model$people$state, model$people$timer, model$people$nodeid,
                                model$people$count, model$nodes$count, Q, R)
  move_count(model$nodes$Q, model$nodes$R, released, model$tick)
}

# ── run: baseline (no intervention) vs quarantine ────────────────────────────────────
base <- run_model(scenario, model = "SIR", nticks = nticks, r0 = r0,
                  infectious_period = inf_duration, seed = 1L)
quar <- run_model(scenario, model = "SIR", nticks = nticks, r0 = r0,
                  infectious_period = inf_duration, seed = 1L, extra_states = "Q",
                  step_update = release_recovered, step_exit = detect_and_isolate)

Ib <- rowSums(base$nodes$I$values()); Iq <- rowSums(quar$nodes$I$values())
cum_b <- cumsum(rowSums(base$nodes$incidence$values()))   # cumulative ever-infected
cum_q <- cumsum(rowSums(quar$nodes$incidence$values()))
cat(sprintf("quarantine: fortnightly testing (every %d d), %.0f%%-sensitive (leaky), %d-day isolation\n",
            test_period, 100 * sensitivity, quarantine_days))
cat(sprintf("  peak infectious:  baseline %s (day %d)  ->  with quarantine %s (day %d)\n",
            format(round(max(Ib)), big.mark = ","), which.max(Ib) - 1L,
            format(round(max(Iq)), big.mark = ","), which.max(Iq) - 1L))
cat(sprintf("  total ever infected: baseline %s  ->  with quarantine %s  (%.0f%% averted)\n",
            format(round(tail(cum_b, 1)), big.mark = ","), format(round(tail(cum_q, 1)), big.mark = ","),
            100 * (1 - tail(cum_q, 1) / tail(cum_b, 1))))
cat(sprintf("  tests run: %d; total detected & isolated: %s\n",
            length(track$day), format(sum(track$detected), big.mark = ",")))

# ── plotting (device-aware: PNG under Rscript, Plots pane in RStudio) ─────────────────
to_png    <- !interactive()
open_png  <- function(path, ...) if (to_png) grDevices::png(path, ...)
close_png <- function() if (to_png) grDevices::dev.off()

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

days <- 0:(nticks - 1L)
Qq   <- rowSums(quar$nodes$Q$values())

open_png(file.path(out_dir, "quarantine.png"), width = 1150, height = 900, res = 120)
op <- graphics::par(mfrow = c(2L, 2L), mar = c(4, 4.5, 2.8, 1))

# (1) Infectious curve, baseline vs quarantine — isolating detected cases flattens it.
matplot(days, cbind(baseline = Ib, quarantine = Iq) / 1e3, type = "l", lty = 1, lwd = 2.5,
        col = c("#d7301f", "#2c7fb8"), xlab = "day", ylab = "infectious (thousands)",
        main = "Infectious: baseline vs quarantine")
graphics::abline(v = test_ticks, col = "grey92")           # testing days
legend("topright", legend = c("baseline", "with quarantine"),
       col = c("#d7301f", "#2c7fb8"), lwd = 2.5, bty = "n")

# (2) Cumulative ever-infected — the gap is cases averted. Incidence is a per-interval flow
# (nticks-1 rows), so plot it against days 1..nticks-1.
matplot(days[-1L], cbind(baseline = cum_b, quarantine = cum_q) / 1e3, type = "l", lty = 1, lwd = 2.5,
        col = c("#d7301f", "#2c7fb8"), xlab = "day", ylab = "cumulative infected (thousands)",
        main = "Cumulative ever-infected (gap = cases averted)")
legend("bottomright", legend = c("baseline", "with quarantine"),
       col = c("#d7301f", "#2c7fb8"), lwd = 2.5, bty = "n")

# (3) The quarantine state Q — pulses up on each fortnightly test day, drains as
# isolated cases are released to R.
plot(days, Qq / 1e3, type = "l", lwd = 2.5, col = "#7570b3",
     xlab = "day", ylab = "in quarantine (thousands)",
     main = "Quarantined (Q): isolated, then released to R")
graphics::abline(v = test_ticks, col = "grey92")

# (4) Detections per test day — the leaky, periodic testing yield.
plot(track$day, track$detected / 1e3, type = "h", lwd = 3, col = "#1b9e77",
     xlab = "day", ylab = "detected per test (thousands)",
     main = sprintf("Cases caught per fortnightly test (%.0f%% sensitive)", 100 * sensitivity))

graphics::par(op)
close_png()
if (to_png) cat(sprintf("wrote quarantine plot to %s\n", file.path(out_dir, "quarantine.png")))
