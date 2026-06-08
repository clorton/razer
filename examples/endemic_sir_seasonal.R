#!/usr/bin/env Rscript

# Endemic SIR with SEASONAL transmission forcing: the two-patch endemic model of
# endemic_sir.R, but with beta multiplied by a gentle annual sinusoid. Seeding near the
# endemic equilibrium leaves R_eff ~ 1, so even a small (+/-10%) wobble produces clear,
# phase-locked ANNUAL epidemic waves.
#
# Built on run_model(): the SIR loop is the runner's; constant-population vital dynamics
# and scheduled importations are added through a `step_exit` callback (with `capacity`
# reserving slots for the imports), and the seasonal multiplier is passed straight to
# run_model's `seasonality` argument (calc_foi multiplies it into the force of infection).
#
# Run from anywhere:  Rscript examples/endemic_sir_seasonal.R

library(razer)

# ── run_endemic_sir: build + run via run_model, report ────────────────────────────────
# As endemic_sir.R, plus a `seasonality` argument (any values_map-broadcastable transmission
# multiplier) forwarded to run_model. Returns (invisibly) the `model` environment.
run_endemic_sir <- function(scenario, network, nticks, inf_duration, r0, cdr, schedule,
                            seasonality = 1, progress = FALSE) {
  if (!inherits(inf_duration, "Distribution"))
    stop("`inf_duration` must be a Distribution (e.g. dist_gamma(2, 4))")

  pops          <- scenario$population
  n             <- sum(pops)
  total_imports <- sum(schedule$count)
  scenario$R    <- pops - round(pops / r0)        # seed S = N/R0, rest immune, I = 0

  sched_tick  <- as.integer(schedule$tick)
  sched_node  <- as.integer(schedule$node)
  sched_count <- as.integer(schedule$count)

  m <- run_model(
    scenario = scenario, model = "SIR", nticks = nticks, r0 = r0,
    infectious_period = inf_duration, network = network, seasonality = seasonality,
    capacity = n + total_imports, seed = 1L, progress = progress,
    init = function(model) {
      nn <- model$nodes$count
      model$nodes$death_rate   <- values_map(cdr / 1000 / 365, nticks, nn)
      model$nodes$births       <- allocate_vector("i32", nticks - 1L, nn)
      model$nodes$deaths       <- allocate_vector("i32", nticks - 1L, nn)
      model$nodes$importations <- allocate_vector("i32", nticks - 1L, nn)
    },
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
scenario <- data.frame(
  name       = c("patch_a", "patch_b"),
  population = c(500000L, 500000L)
)
n_nodes <- nrow(scenario)
network <- matrix(c(0.00, 0.01,
                    0.01, 0.00), nrow = n_nodes, byrow = TRUE)

r0           <- 3
inf_duration <- dist_gamma(2, 4)     # infectious period (mean shape*scale = 8 ticks)
cdr          <- 20                   # crude death rate (annual per 1,000)
nticks       <- 3650L                # ten years of daily steps

# ── seasonal transmission forcing ─────────────────────────────────────────────────
# A gentle annual sinusoid on the transmission MULTIPLIER: beta is scaled by
# 1 + amplitude*sin(2*pi*day/365). We use +/-10% — the smallest forcing that yields clear,
# phase-locked annual waves whose off-season trough never collapses to zero (larger
# amplitudes occasionally faded out a patch). `seq_len(nticks) - 1L` is the 0-based day
# index; run_model's values_map broadcasts this length-nticks vector across both patches.
seasonal_amplitude <- 0.10
seasonality        <- 1 + seasonal_amplitude * sin(2 * pi * (seq_len(nticks) - 1L) / 365)

# Importation schedule: 10 new infectious cases into EACH patch every 30 ticks — a small
# endemic floor so an off-season trough can't permanently fade out.
import_ticks   <- seq(0L, nticks - 2L, by = 30L)
schedule       <- expand.grid(tick = import_ticks, node = seq_len(n_nodes) - 1L)
schedule$count <- 10L

# ── run ───────────────────────────────────────────────────────────────────────────
timing <- system.time(
  result <- run_endemic_sir(scenario, network, nticks, inf_duration, r0, cdr, schedule,
                            seasonality = seasonality, progress = TRUE))
cat(sprintf("run_endemic_sir (seasonal) completed in %.3f s\n", timing[["elapsed"]]))

# ── plot the seasonal forcing and the annual epidemic response ──────────────────────
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

I_total <- rowSums(result$nodes$I$values())
years   <- (seq_len(nticks) - 1L) / 365

plot_path <- file.path(out_dir, "endemic_sir_seasonal.png")
open_png(plot_path, width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))
plot(years, seasonality, type = "l", col = "darkorange", lwd = 2,
     xlab = "time (years)", ylab = "beta multiplier",
     main = sprintf("Seasonal transmission forcing (+/-%.0f%%)", 100 * seasonal_amplitude))
graphics::abline(v = 0:10, col = "grey85"); graphics::abline(h = 1, col = "grey60", lty = 3)
plot(years, I_total, type = "l", col = "firebrick", lwd = 2,
     xlab = "time (years)", ylab = "agents",
     main = "Infectious prevalence (summed over patches): annual epidemic waves")
graphics::abline(v = 0:10, col = "grey85")
graphics::par(op)
close_png()
if (to_png) cat(sprintf("wrote seasonal trajectory plot to %s\n", plot_path))

# Per-year peak and trough of total prevalence: the run should oscillate annually without
# the trough collapsing to zero. `tapply(x, g, f)` applies `f` to `x` split by group `g`.
yr      <- (seq_along(I_total) - 1L) %/% 365L
peaks   <- tapply(I_total, yr, max)
troughs <- tapply(I_total, yr, min)
cat("per-year infectious peak:   ", paste(as.integer(peaks),   collapse = ", "), "\n", sep = "")
cat("per-year infectious trough: ", paste(as.integer(troughs), collapse = ", "), "\n", sep = "")
