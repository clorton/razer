#!/usr/bin/env Rscript

# Quarantine / case isolation, via a user-defined "Q" agent state. Each tick an intervention
# "tests" the infectious population and, with probability `detection_rate`, moves a detected
# infectious agent I -> Q. A quarantined agent is isolated: Q is not the I state, and the
# force of infection reads only the I census, so quarantined cases NO LONGER TRANSMIT. After
# a fixed `quarantine_days` isolation period a step_update callback converts Q -> R (the
# agent has recovered) using the generic step_timer_expire kernel.
#
# Q is registered with run_model(extra_states = "Q"); the disease kernels leave Q agents
# untouched, so both the isolation (I->Q) and the release (Q->R) are composed from callbacks
# with no kernel change. We run the same epidemic twice — without and with quarantine — to
# show the intervention flattening the curve.
#
# Run from anywhere:  Rscript examples/quarantine.R

library(razer)

set.seed(20260608L)

# ── parameters ────────────────────────────────────────────────────────────────────
N               <- 500000L
nticks          <- 200L
r0              <- 2.5
inf_duration    <- dist_gamma(2, 4)     # infectious period (mean 8 ticks)
detection_rate  <- 0.08                  # per-tick probability an infectious agent is detected
quarantine_days <- 10L                   # isolation period, then Q -> R
scenario        <- data.frame(population = N, I = 100L)

# ── intervention callbacks ──────────────────────────────────────────────────────────
# detect_and_isolate (step_exit): each tick, flip each infectious agent to Q with prob
# `detection_rate`, set its isolation timer, and apply the per-node I -> Q census delta.
detect_and_isolate <- function(model) {
  I <- model$states[["I"]]; Q <- model$states[["Q"]]
  cnt   <- model$people$count
  state <- model$people$state$values()
  node  <- model$people$nodeid$values()
  infectious <- which(state[seq_len(cnt)] == I)
  detected   <- infectious[stats::runif(length(infectious)) < detection_rate]
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
cat(sprintf("quarantine: detection %.0f%%/day, %d-day isolation\n", 100 * detection_rate, quarantine_days))
cat(sprintf("  peak infectious:  baseline %s (day %d)  ->  with quarantine %s (day %d)\n",
            format(round(max(Ib)), big.mark = ","), which.max(Ib) - 1L,
            format(round(max(Iq)), big.mark = ","), which.max(Iq) - 1L))
cat(sprintf("  ever infected (final R + Q): baseline %s  ->  with quarantine %s\n",
            format(base$nodes$R$values()[nticks, 1L], big.mark = ","),
            format(quar$nodes$R$values()[nticks, 1L] + quar$nodes$Q$values()[nticks, 1L], big.mark = ",")))

# ── plot ─────────────────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

days <- 0:(nticks - 1L)
Qq   <- rowSums(quar$nodes$Q$values())

grDevices::png(file.path(out_dir, "quarantine.png"), width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))

# Top: infectious curve, baseline vs quarantine — isolating detected cases flattens it.
matplot(days, cbind(baseline = Ib, quarantine = Iq) / 1e3, type = "l", lty = 1, lwd = 2.5,
        col = c("#d7301f", "#2c7fb8"), xlab = "day", ylab = "infectious (thousands)",
        main = sprintf("Infectious curve: baseline vs quarantine (%.0f%%/day detection, %d-day isolation)",
                       100 * detection_rate, quarantine_days))
legend("topright", legend = c("baseline (no quarantine)", "with quarantine"),
       col = c("#d7301f", "#2c7fb8"), lwd = 2.5, bty = "n")

# Bottom: the quarantine compartment Q — infectious agents currently isolated.
plot(days, Qq / 1e3, type = "l", lwd = 2.5, col = "#7570b3",
     xlab = "day", ylab = "in quarantine (thousands)",
     main = "Quarantined (Q): detected infectious agents, isolated then released to R")
grDevices::dev.off()
cat(sprintf("wrote quarantine plot to %s\n", file.path(out_dir, "quarantine.png")))
