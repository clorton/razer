#!/usr/bin/env Rscript

# Supplemental Immunization Activities (SIA) — periodic vaccination campaigns. A campaign
# runs on a configurable schedule (here once a year) and, in the nodes it targets,
# vaccinates each susceptible with probability `coverage`, moving them S -> V. V is a
# user-defined agent state registered with run_model(extra_states = "V"); the disease
# kernels leave V agents untouched (protected — never infected, never transitioned), and
# run_model tracks the V census for us. There is NO waning here, so V is permanent and the
# campaign is the only thing that creates it — implemented purely as a step_exit callback
# (no step_update). See sia_campaigns_waning.R for the waning variant.
#
# To make the node-targeting visible, the three patches are INDEPENDENT (no coupling) and
# the campaign targets only nodes 0 and 1; node 2 is the unvaccinated control.
#
# Run from anywhere:  Rscript examples/sia_campaigns.R

library(razer)

set.seed(20260608L)   # R-side campaign draws (kernel RNG is seeded separately, below)

# ── parameters ────────────────────────────────────────────────────────────────────
pop_per_node <- 300000L
n_nodes      <- 3L
years        <- 4L
nticks       <- years * 365L
r0           <- 2
inf_duration <- dist_gamma(2, 4)              # infectious period (mean 8 ticks)

# Campaign schedule + targeting (the configurable knobs):
campaign_period <- 365L                       # once a year
campaign_start  <- 30L                        # first campaign early (preventive, before the wave)
campaign_ticks  <- seq(campaign_start, nticks - 2L, by = campaign_period)
affected_nodes  <- c(0L, 1L)                  # 0-based node ids the campaign reaches (not node 2)
coverage        <- 0.6                         # per-susceptible vaccination probability

scenario <- data.frame(population = rep(pop_per_node, n_nodes), I = rep(50L, n_nodes))
network  <- matrix(0, n_nodes, n_nodes)        # independent patches: targeting is crisp

# ── campaign callback: probabilistic S -> V in the affected nodes, on schedule ────────
# Runs at the end of each tick (step_exit); model$tick is the 0-based interval index. On a
# campaign tick it flips each eligible susceptible to V with probability `coverage`, then
# applies the per-node S -> V census delta with move_count.
sia_campaign <- function(model) {
  if (!(model$tick %in% campaign_ticks)) return(invisible(NULL))
  V <- model$states[["V"]]; S <- model$states[["S"]]
  cnt   <- model$people$count
  state <- model$people$state$values()
  node  <- model$people$nodeid$values()
  eligible <- which(state[seq_len(cnt)] == S & node[seq_len(cnt)] %in% affected_nodes)
  vacc     <- eligible[stats::runif(length(eligible)) < coverage]   # probabilistic coverage
  if (!length(vacc)) return(invisible(NULL))
  state[vacc] <- V
  model$people$state$set(state)
  counts <- tabulate(node[vacc] + 1L, model$nodes$count)            # vaccinated per node
  move_count(model$nodes$S, model$nodes$V, counts, model$tick)      # S -> V at slice tick+1
}

# ── run ─────────────────────────────────────────────────────────────────────────────
m <- run_model(scenario, model = "SIR", nticks = nticks, r0 = r0,
               infectious_period = inf_duration, network = network, seed = 1L,
               extra_states = "V", step_exit = sia_campaign)

cat(sprintf("SIA campaigns: %d patches x %s, %d years, %d annual campaigns at %.0f%% coverage in nodes {%s}\n",
            n_nodes, format(pop_per_node, big.mark = ","), years, length(campaign_ticks),
            100 * coverage, paste(affected_nodes, collapse = ", ")))
Vend <- m$nodes$V$values()[nticks, ]
for (k in seq_len(n_nodes))
  cat(sprintf("  node %d: final vaccinated (V) = %s,  final recovered (R) = %s\n",
              k - 1L, format(Vend[k], big.mark = ","),
              format(m$nodes$R$values()[nticks, k], big.mark = ",")))

# ── plots ─────────────────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

years_axis <- (seq_len(nticks) - 1L) / 365
node_col   <- c("#1b9e77", "#7570b3", "#d95f02")     # nodes 0, 1, 2
Vm <- m$nodes$V$values() / 1e3                       # thousands
Im <- m$nodes$I$values() / 1e3

grDevices::png(file.path(out_dir, "sia_campaigns.png"), width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))

# Top: vaccinated (V) per node — the two targeted nodes step up at each annual campaign,
# the untargeted node 2 stays flat at zero.
matplot(years_axis, Vm, type = "l", lty = 1, lwd = 2, col = node_col,
        xlab = "year", ylab = "vaccinated (thousands)",
        main = sprintf("Vaccinated (V) per node — annual SIAs at %.0f%% coverage in nodes 0 & 1", 100 * coverage))
graphics::abline(v = campaign_ticks / 365, col = "grey85", lty = 3)
legend("topleft", legend = c("node 0 (targeted)", "node 1 (targeted)", "node 2 (control)"),
       col = node_col, lwd = 2, bty = "n")

# Bottom: infectious (I) per node — the vaccinated nodes are protected, node 2 burns through.
matplot(years_axis, Im, type = "l", lty = 1, lwd = 2, col = node_col,
        xlab = "year", ylab = "infectious (thousands)",
        main = "Infectious (I) per node — the unvaccinated control sustains more transmission")
graphics::abline(v = campaign_ticks / 365, col = "grey85", lty = 3)
legend("topright", legend = c("node 0", "node 1", "node 2 (control)"),
       col = node_col, lwd = 2, bty = "n")
graphics::par(op)
grDevices::dev.off()
cat(sprintf("wrote SIA plot to %s\n", file.path(out_dir, "sia_campaigns.png")))
