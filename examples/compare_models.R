#!/usr/bin/env Rscript

# Compare four models from the menagerie — SIR, SIRS, SEIR, SEIRS — on the SAME population,
# duration, and transmission parameters, so the only differences are the state
# structure (an exposed/latent E stage) and waning immunity (R -> S).
#
#   * population     1,000,000 (single well-mixed node), 100 initially infectious
#   * duration       365 ticks (days)
#   * R0 / beta       2.5 / ~0.35  (run_model derives beta = R0 / mean(infectious_period))
#   * infectious      gamma(shape = 140, scale = 0.05)  -> mean 7 days (so beta*D ~= 2.5)
#   * incubation      gamma(shape = 140, scale = 0.05)  -> mean 7 days (SEIR/SEIRS only)
#   * immunity        60 days (SIRS/SEIRS only; not part of the shared comparison, chosen
#                     so waning resupplies susceptibles within the year)
#
# Each trajectory is colored by STATE (S blue, E orange, I red, R green) and styled by
# MODEL (line type), so e.g. all four S curves are blue but distinguishable by line style.
#
# Run from anywhere:  Rscript examples/compare_models.R
# Output PNGs are written next to this script in examples/output/.

library(razer)

set.seed(20260608L)

# ── shared parameters ─────────────────────────────────────────────────────────────
N        <- 1000000L
nticks   <- 365L
R0       <- 2.5
inf_dist <- dist_gamma(140, 0.05)        # infectious period: mean shape*scale = 7 days
inc_dist <- dist_gamma(140, 0.05)        # incubation period (E stage): same mean 7 days
imm_days <- 60                           # waning-immunity period (SIRS / SEIRS only)
scenario <- data.frame(population = N, I = 100L)
models   <- c("SIR", "SIRS", "SEIR", "SEIRS")

# ── run each model with exactly the parameters it needs ─────────────────────────────
# Passing only the relevant durations avoids run_model's "ignored parameter" warnings:
# E models get an incubation period; waning models (SIRS/SEIRS) get an immunity period.
run_one <- function(model) {
  args <- list(scenario = scenario, model = model, nticks = nticks, r0 = R0,
               infectious_period = inf_dist, seed = 1L)
  if (grepl("E", model))             args$incubation_period <- inc_dist
  if (model %in% c("SIRS", "SEIRS")) args$immunity_period   <- imm_days
  do.call(run_model, args)
}
results <- lapply(models, run_one)
names(results) <- models

# Pull each present state's trajectory (single node -> column 1 of the census matrix).
traj <- lapply(results, function(m) {
  comp <- list(S = m$nodes$S$values()[, 1], I = m$nodes$I$values()[, 1])
  if (!is.null(m$nodes$E)) comp$E <- m$nodes$E$values()[, 1]
  if (!is.null(m$nodes$R)) comp$R <- m$nodes$R$values()[, 1]
  comp
})

# ── styling: colour = state, line type = model ───────────────────────────────
col_comp <- c(S = "#1f78b4", E = "#ff7f00", I = "#e31a1c", R = "#33a02c")   # blue/orange/red/green
lty_model <- c(SIR = 1, SIRS = 2, SEIR = 3, SEIRS = 5)                      # solid/dashed/dotted/longdash
ticks <- 0:(nticks - 1L)

# ── output location (next to this script) ───────────────────────────────────────────
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

# ── plot 1: all trajectories overlaid (colour = state, line type = model) ─────
open_png(file.path(out_dir, "compare_models.png"), width = 1150, height = 760, res = 120)
graphics::par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = c(0, nticks - 1L), ylim = c(0, N / 1e6),
     xlab = "day", ylab = "agents (millions)",
     main = sprintf("SIR / SIRS / SEIR / SEIRS  (N = %s, R0 = %.1f, D = 7 d)",
                    format(N, big.mark = ","), R0))
for (model in models)
  for (comp in names(traj[[model]]))
    graphics::lines(ticks, traj[[model]][[comp]] / 1e6,
                    col = col_comp[comp], lty = lty_model[model], lwd = 2)
legend("top", title = "state", legend = names(col_comp), col = col_comp,
       lwd = 2, lty = 1, bty = "n", horiz = TRUE)
legend("right", title = "model", legend = names(lty_model), col = "grey30",
       lwd = 2, lty = lty_model, bty = "n")
close_png()

# ── plot 2: one panel per state (the four models compared within each) ────────
# Clearer than the overlay when curves cross; same colour/line-type scheme.
open_png(file.path(out_dir, "compare_models_panels.png"), width = 1150, height = 900, res = 120)
op <- graphics::par(mfrow = c(2L, 2L), mar = c(4, 4.2, 2.5, 1))
for (comp in c("S", "E", "I", "R")) {
  have <- models[vapply(models, function(m) !is.null(traj[[m]][[comp]]), logical(1))]
  ymax <- max(vapply(have, function(m) max(traj[[m]][[comp]]), numeric(1))) / 1e6
  plot(NA, xlim = c(0, nticks - 1L), ylim = c(0, ymax),
       xlab = "day", ylab = "agents (millions)",
       main = sprintf("%s state", comp), col.main = col_comp[comp])
  for (model in have)
    graphics::lines(ticks, traj[[model]][[comp]] / 1e6,
                    col = col_comp[comp], lty = lty_model[model], lwd = 2)
  legend("topright", have, col = "grey30", lwd = 2, lty = lty_model[have], bty = "n")
}
graphics::par(op)
close_png()

# ── summary ─────────────────────────────────────────────────────────────────────────
cat("model    peak I (day)     final S      final R\n")
for (model in models) {
  I <- traj[[model]]$I; S <- traj[[model]]$S; R <- traj[[model]]$R
  cat(sprintf("%-7s  %s (%3d)   %s   %s\n", model,
              format(round(max(I)), big.mark = ",", width = 9), which.max(I) - 1L,
              format(round(tail(S, 1)), big.mark = ",", width = 9),
              if (is.null(R)) "       --" else format(round(tail(R, 1)), big.mark = ",", width = 9)))
}
if (to_png) cat(sprintf("\nPlots written to: %s\n", normalizePath(out_dir)))
