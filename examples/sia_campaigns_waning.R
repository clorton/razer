#!/usr/bin/env Rscript

# Supplemental Immunization Activities (SIA) WITH waning vaccine-derived immunity. Same
# periodic, node-targeted, probabilistic campaign as sia_campaigns.R (S -> V at `coverage`
# in the affected nodes once a year), but now protection is TEMPORARY: each vaccination
# sets a per-agent waning timer, and a step_update callback runs the generic
# step_timer_expire(V -> S) kernel every tick to return agents to susceptible when their
# immunity lapses. Nothing in the disease kernels changes — V waning is composed entirely
# from the existing kernel + the model$states code map.
#
# The visible difference from the no-waning example: the vaccinated (V) count SAWTOOTHS —
# it jumps at each annual campaign and decays between — instead of staircasing up forever.
# Periodic campaigns are needed precisely because immunity wanes.
#
# Run from anywhere:  Rscript examples/sia_campaigns_waning.R

library(razer)

set.seed(20260608L)

# ── parameters ────────────────────────────────────────────────────────────────────
pop_per_node <- 300000L
n_nodes      <- 3L
years        <- 6L
nticks       <- years * 365L
r0           <- 2
inf_duration <- dist_gamma(2, 4)              # infectious period (mean 8 ticks)

campaign_period <- 365L                       # once a year
campaign_start  <- 30L
campaign_ticks  <- seq(campaign_start, nticks - 2L, by = campaign_period)
affected_nodes  <- c(0L, 1L)                  # node 2 is the unvaccinated control
coverage        <- 0.6

# NEW: vaccine-derived immunity wanes after ~2 years (per-agent draw, so protection lapses
# gradually rather than at a sharp cliff). This duration is the single new knob vs the
# no-waning example.
waning_period <- dist_normal(2 * 365, 90)     # mean ~730 days, sd 90

scenario <- data.frame(population = rep(pop_per_node, n_nodes), I = rep(50L, n_nodes))
network  <- matrix(0, n_nodes, n_nodes)

# ── campaign callback: S -> V at `coverage`, SETTING a waning timer on each vaccinee ──
sia_campaign <- function(model) {
  if (!(model$tick %in% campaign_ticks)) return(invisible(NULL))
  V <- model$states[["V"]]; S <- model$states[["S"]]
  cnt   <- model$people$count
  state <- model$people$state$values()
  node  <- model$people$nodeid$values()
  eligible <- which(state[seq_len(cnt)] == S & node[seq_len(cnt)] %in% affected_nodes)
  vacc     <- eligible[stats::runif(length(eligible)) < coverage]
  if (!length(vacc)) return(invisible(NULL))
  state[vacc] <- V
  timer <- model$people$timer$values()
  timer[vacc] <- pmax(1L, pmin(65535L, as.integer(round(waning_period$sample_n(length(vacc))))))
  model$people$state$set(state)
  model$people$timer$set(timer)               # the V -> S countdown
  counts <- tabulate(node[vacc] + 1L, model$nodes$count)
  move_count(model$nodes$S, model$nodes$V, counts, model$tick)
}

# ── waning callback: V -> S on timer expiry, via the generic step_timer_expire kernel ─
# This is the only addition over sia_campaigns.R. step_timer_expire decrements every V
# agent's timer and, on reaching 0, moves it back to S (returning per-node counts); no new
# kernel and no change to the disease step. We tally the waning events for the summary.
track <- new.env(); track$waned <- 0
vaccine_waning <- function(model) {
  V <- model$states[["V"]]; S <- model$states[["S"]]
  waned <- step_timer_expire(model$people$state, model$people$timer, model$people$nodeid,
                             model$people$count, model$nodes$count, V, S)
  move_count(model$nodes$V, model$nodes$S, waned, model$tick)
  track$waned <- track$waned + sum(waned)
}

# ── run ─────────────────────────────────────────────────────────────────────────────
m <- run_model(scenario, model = "SIR", nticks = nticks, r0 = r0,
               infectious_period = inf_duration, network = network, seed = 1L,
               extra_states = "V", step_update = vaccine_waning, step_exit = sia_campaign)

cat(sprintf("SIA + waning: %d patches, %d years, annual campaigns (%.0f%% coverage, nodes {%s}), mean waning ~%.0f d\n",
            n_nodes, years, 100 * coverage, paste(affected_nodes, collapse = ", "), 2 * 365))
cat(sprintf("  total V->S waning events over the run: %s; final vaccinated (V) per node: %s\n",
            format(round(track$waned), big.mark = ","),
            paste(format(m$nodes$V$values()[nticks, ], big.mark = ","), collapse = ", ")))

# ── plots ─────────────────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

years_axis <- (seq_len(nticks) - 1L) / 365
node_col   <- c("#1b9e77", "#7570b3", "#d95f02")
Vm <- m$nodes$V$values() / 1e3
Sm <- m$nodes$S$values() / 1e3

grDevices::png(file.path(out_dir, "sia_campaigns_waning.png"), width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))

# Top: vaccinated (V) per node — a SAWTOOTH (campaign jumps, then waning decay) in the two
# targeted nodes, vs the permanent staircase of the no-waning example.
matplot(years_axis, Vm, type = "l", lty = 1, lwd = 2, col = node_col,
        xlab = "year", ylab = "vaccinated (thousands)",
        main = "Vaccinated (V) per node — annual SIAs with ~2-year waning (sawtooth, not staircase)")
graphics::abline(v = campaign_ticks / 365, col = "grey85", lty = 3)
legend("topleft", legend = c("node 0 (targeted)", "node 1 (targeted)", "node 2 (control)"),
       col = node_col, lwd = 2, bty = "n")

# Bottom: susceptibles (S) per node — they dip at each campaign and refill as V wanes back,
# so the campaign must be repeated to hold coverage.
matplot(years_axis, Sm, type = "l", lty = 1, lwd = 2, col = node_col,
        xlab = "year", ylab = "susceptible (thousands)",
        main = "Susceptible (S) per node — refilled by waning between campaigns")
graphics::abline(v = campaign_ticks / 365, col = "grey85", lty = 3)
legend("right", legend = c("node 0", "node 1", "node 2 (control)"),
       col = node_col, lwd = 2, bty = "n")
graphics::par(op)
grDevices::dev.off()
cat(sprintf("wrote SIA-with-waning plot to %s\n", file.path(out_dir, "sia_campaigns_waning.png")))
