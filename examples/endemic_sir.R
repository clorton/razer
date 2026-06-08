#!/usr/bin/env Rscript

# Endemic SIR example: a two-patch metapopulation driven to ENDEMICITY by constant-
# population vital dynamics (turnover continually resupplies susceptibles) plus periodic
# IMPORTATIONS that re-spark transmission so a stochastic fade-out doesn't end the epidemic.
#
# Built on the high-level runner run_model(), which owns the disease loop. Two things go
# beyond the closed-population menagerie and are added via run_model's hooks:
#   * `capacity` — run_model is given room ABOVE the initial population so import_infections
#     has reserved agent slots to activate (a closed run, capacity == population, could not
#     grow). Sized as initial population + total scheduled imports.
#   * a `step_exit` callback — constant_pop_vitals_sir (births = deaths, reborn susceptible)
#     then import_infections (activate reserved slots as new infectious cases on schedule),
#     run after transmission each tick, exactly as the hand-wired loop did.
#
# Run from anywhere:  Rscript examples/endemic_sir.R

library(razer)

# ── run_endemic_sir: build + run via run_model, report ────────────────────────────────
# `scenario` has one row per patch (name, population); `network` the N x N coupling matrix;
# `nticks` recorded daily states; `inf_duration` the infectious-period Distribution; `r0`
# the basic reproduction number; `cdr` the crude death rate (annual per 1,000); `schedule`
# a data.frame of importations (integer 0-based `tick` / `node`, integer `count`). Returns
# (invisibly) the `model` environment.
run_endemic_sir <- function(scenario, network, nticks, inf_duration, r0, cdr, schedule,
                            progress = FALSE) {
  if (!inherits(inf_duration, "Distribution"))
    stop("`inf_duration` must be a Distribution (e.g. dist_gamma(2, 4))")

  pops          <- scenario$population
  n             <- sum(pops)
  total_imports <- sum(schedule$count)

  # Seed at the endemic equilibrium: susceptible fraction 1/R0, the rest immune, I = 0
  # (sparked by imports). run_model seeds the per-node `R` column as recovered and makes S
  # the remainder, so set R = population - round(population / R0).
  scenario$R <- pops - round(pops / r0)

  # Importation schedule as three integer vectors the Rust kernel scans per tick; captured
  # by the step_exit closure below.
  sched_tick  <- as.integer(schedule$tick)
  sched_node  <- as.integer(schedule$node)
  sched_count <- as.integer(schedule$count)

  # CAPACITY above the initial population leaves reserved slots [n, n + total_imports) that
  # import_infections activates. Vital dynamics are constant-population (no growth), so only
  # the imports need the headroom.
  m <- run_model(
    scenario = scenario, model = "SIR", nticks = nticks, r0 = r0,
    infectious_period = inf_duration, network = network,
    capacity = n + total_imports, seed = 1L, progress = progress,
    init = function(model) {
      nn <- model$nodes$count
      model$nodes$death_rate   <- values_map(cdr / 1000 / 365, nticks, nn)  # daily death hazard
      model$nodes$births       <- allocate_vector("i32", nticks - 1L, nn)
      model$nodes$deaths       <- allocate_vector("i32", nticks - 1L, nn)
      model$nodes$importations <- allocate_vector("i32", nticks - 1L, nn)
    },
    # After transmission each tick: constant-pop births/deaths, then imports (last, so the
    # imported cases seed the NEXT tick's force of infection and the next carry-forward
    # folds them into N). import_infections returns the grown live count.
    step_exit = function(model) {
      constant_pop_vitals_sir(model$people$state, model$people$timer, model$people$nodeid,
        model$people$count, model$nodes$death_rate, model$nodes$S, model$nodes$I, model$nodes$R,
        model$nodes$births, model$nodes$deaths, model$tick)
      model$people$count <- import_infections(
        model$people$state, model$people$timer, model$people$nodeid, model$people$count,
        model$nodes$I, model$nodes$importations, sched_tick, sched_node, sched_count,
        inf_duration, model$tick)
    })

  final <- nticks
  cat(sprintf(paste0("run_endemic_sir: %d patches, %d -> %d agents (%d imported), %d ticks; ",
                     "final S=%d I=%d R=%d N=%d; incidence=%d, recoveries=%d, births=deaths=%d.\n"),
              m$nodes$count, n, m$people$count, total_imports, nticks,
              sum(m$nodes$S$values()[final, ]), sum(m$nodes$I$values()[final, ]),
              sum(m$nodes$R$values()[final, ]), sum(m$nodes$N$values()[final, ]),
              sum(m$nodes$incidence$values()), sum(m$nodes$recovery$values()),
              sum(m$nodes$births$values())))
  invisible(m)
}

# ── setup ───────────────────────────────────────────────────────────────────────────
# Two equal patches of 500,000 agents each.
scenario <- data.frame(
  name       = c("patch_a", "patch_b"),
  population = c(500000L, 500000L)
)
n_nodes <- nrow(scenario)

# A small symmetric cross-coupling so the patches exchange 1% of their force of infection.
network <- matrix(c(0.00, 0.01,
                    0.01, 0.00), nrow = n_nodes, byrow = TRUE)

r0           <- 3
inf_duration <- dist_gamma(2, 4)     # infectious period (mean shape*scale = 8 ticks)
cdr          <- 20                   # crude death rate (annual per 1,000): turnover sustains endemicity
nticks       <- 3650L                # ten years of daily steps

# Importation schedule: 10 new infectious cases into EACH patch every 30 ticks.
import_ticks   <- seq(0L, nticks - 2L, by = 30L)
schedule       <- expand.grid(tick = import_ticks, node = seq_len(n_nodes) - 1L)
schedule$count <- 10L

# ── run ───────────────────────────────────────────────────────────────────────────
timing <- system.time(
  result <- run_endemic_sir(scenario, network, nticks, inf_duration, r0, cdr, schedule,
                            progress = TRUE))
cat(sprintf("run_endemic_sir completed in %.3f s\n", timing[["elapsed"]]))

# ── plot the S, I, R channels over time ─────────────────────────────────────────────
# `$values()` returns each census buffer as an nticks x n_nodes matrix; `rowSums` totals
# over patches. Time is in years (tick / 365).
S_total <- rowSums(result$nodes$S$values())
I_total <- rowSums(result$nodes$I$values())
R_total <- rowSums(result$nodes$R$values())
years   <- (seq_len(nticks) - 1L) / 365

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# S and R are ~10^5 while I is ~10^1, so a shared linear axis would flatten I — stack two
# panels and give I its own.
plot_path <- file.path(out_dir, "endemic_sir_SIR.png")
grDevices::png(plot_path, width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))
plot(years, S_total, type = "l", col = "steelblue", lwd = 2,
     ylim = range(0, S_total, R_total),
     xlab = "time (years)", ylab = "agents",
     main = "Endemic SIR: susceptible & recovered (summed over patches)")
lines(years, R_total, col = "darkgreen", lwd = 2)
graphics::abline(h = sum(scenario$population) / r0, col = "steelblue", lty = 3)  # S* = N/R0
legend("right", legend = c("S", "R", "N / R0"), bty = "n",
       col = c("steelblue", "darkgreen", "steelblue"),
       lwd = c(2, 2, 1), lty = c(1, 1, 3))
plot(years, I_total, type = "l", col = "firebrick", lwd = 2,
     xlab = "time (years)", ylab = "agents",
     main = "Endemic SIR: infectious prevalence")
legend("topright", legend = "I", bty = "n", col = "firebrick", lwd = 2)
graphics::par(op)
grDevices::dev.off()
cat(sprintf("wrote S/I/R trajectory plot to %s\n", plot_path))
