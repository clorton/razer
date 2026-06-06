#!/usr/bin/env Rscript

# Seasonally-forced endemic SIR example. This is endemic_sir.R (a two-patch
# metapopulation kept endemic by constant-population vital turnover plus periodic
# importations) with ONE addition: the transmission coefficient is modulated by a
# gentle annual sinusoid. The interesting result is how SMALL a forcing it takes to
# produce pronounced ANNUAL epidemic waves: because the population is seeded at the
# endemic susceptible fraction S ≈ N/R0, the effective reproduction number sits right
# at ~1 (critically poised), so even a +/-10% seasonal wobble in beta is amplified into
# large, phase-locked yearly outbreaks (peaks recur around the same time each year).
#
# The three pieces beyond simple_sir.R are the same as endemic_sir.R:
#   * import_infections()      — activates RESERVED agent slots as new infectious cases;
#   * carry_forward_states()   — carries S/I/R forward and totals them into N;
#   * constant_pop_vitals_sir()— births = deaths, reborn susceptible, keeping S+I+R=N.
# The only new ingredient is a per-tick `seasonality` grid passed to calc_foi().
#
# Run from anywhere:  Rscript examples/endemic_sir_seasonal.R

# `library(pkg)` attaches a package so its exported names resolve unqualified (like a
# wildcard `import`); the package name is given bare, not as a string.
library(razer)

# ── run_endemic_sir: build the model, run it, report ────────────────────────────────
# Identical to endemic_sir.R's runner except for the `seasonality` argument: a
# transmission MULTIPLIER accepted in any broadcastable form (a scalar, a per-node
# vector, a per-tick vector of length `nticks`, or a full nticks x nnodes matrix);
# `values_map` expands it to the grid calc_foi() multiplies into the local force of
# infection. `scenario` is a data.frame (name, population) per patch; `network` the
# N x N spatial coupling; `nticks` the recorded daily states; `inf_duration` the
# infectious-period Distribution; `r0` the basic reproduction number (beta = r0/mean);
# `cdr` the crude death rate (annual deaths per 1000); `schedule` the importation table
# (integer columns tick, node, count, all 0-based for tick/node); `progress` draws a
# text progress bar. Returns (invisibly) a list with the `people` and `nodes` envs.
run_endemic_sir <- function(scenario, network, nticks, inf_duration, r0, cdr, schedule,
                            seasonality = 1, progress = FALSE) {
  if (!inherits(inf_duration, "Distribution"))
    stop("`inf_duration` must be a Distribution (e.g. dist_gamma(2, 4))")

  # beta from R0: R0 = beta * mean infectious duration, so beta = R0 / mean. Estimate
  # the mean by sampling the distribution (a large batch; `mean()` averages it).
  beta    <- r0 / mean(inf_duration$sample_n(100000L))
  pops    <- scenario$population
  n       <- sum(pops)              # initial live-agent count (sum over patches)
  n_nodes <- nrow(scenario)

  # CAPACITY above COUNT: every scheduled importation activates a fresh agent slot, so
  # the per-agent arrays must be allocated with room for all of them up front. `count`
  # is the number of currently-live agents; `capacity` the allocated array length.
  total_imports <- sum(schedule$count)
  capacity      <- n + total_imports

  # `new.env()` makes a fresh ENVIRONMENT: unlike ordinary R values (copy-on-modify),
  # an environment has reference semantics — assigning into it mutates the same object
  # everyone holds (a mutable struct). `env$name <- value` assigns a member by name.
  people <- new.env()
  people$count    <- n
  people$capacity <- capacity
  # `allocate_scalar()` (Rust) returns an opaque Column handle over a Rust-owned,
  # zero-filled buffer of the given dtype/length. Sized to CAPACITY so reserved slots
  # [count, capacity) exist for import_infections() to fill. state defaults to 0 = S.
  people$state    <- allocate_scalar("u8",  capacity)
  people$nodeid   <- allocate_scalar("u16", capacity)
  people$timer    <- allocate_scalar("u8",  capacity)

  # nodeid is 0-BASED (0..N-1), matching the Rust kernels (they index per-node arrays
  # directly by it). `rep(ids, times = pops)` repeats each node id by its patch's
  # population; reserved slots past `n` stay 0 until an import sets them. `$set()`
  # writes the WHOLE capacity-length buffer in one copy, so pad to `capacity`.
  ids <- seq_len(n_nodes) - 1L
  people$nodeid$set(c(rep(ids, times = pops), rep(0L, total_imports)))

  # ── per-patch record (census buffers + flow reports + driver grids) ───────────────
  # The S/I/R CENSUS buffers are nticks x n_nodes (state 0 the initial condition, state
  # nticks-1 the final; maintained incrementally — each tick carries column t forward to
  # t+1 and applies only deltas). N (population / FOI denominator) is the same shape,
  # recomputed each tick as S+I+R so it tracks the agents added by imports. The
  # per-interval FLOW reports are (nticks-1) x n_nodes. All are time-major (each tick's
  # per-node row is contiguous).
  nodes <- new.env()
  nodes$count        <- n_nodes
  nodes$S            <- allocate_vector("i32", nticks,      n_nodes)
  nodes$I            <- allocate_vector("i32", nticks,      n_nodes)
  nodes$R            <- allocate_vector("i32", nticks,      n_nodes)
  nodes$N            <- allocate_vector("i32", nticks,      n_nodes)
  nodes$foi          <- allocate_vector("f64", nticks - 1L, n_nodes)
  nodes$recoveries   <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$incidence    <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$births       <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$deaths       <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$importations <- allocate_vector("i32", nticks - 1L, n_nodes)
  # Transmission / vital-dynamics driver grids (nticks x n_nodes f64) built by
  # values_map. `seasonality` is now a per-tick sinusoid rather than the constant 1 of
  # endemic_sir.R. The death rate is the per-node daily death HAZARD (CDR per 1000 per
  # year -> daily).
  nodes$beta         <- values_map(beta,             nticks, n_nodes)
  nodes$seasonality  <- values_map(seasonality,      nticks, n_nodes)
  nodes$death_rate   <- values_map(cdr / 1000 / 365, nticks, n_nodes)

  # ── seed the initial condition: 1/R0 susceptible, the rest recovered, I = 0 ───────
  # Seeding at the endemic susceptible fraction 1/R0 is exactly what makes this model
  # so sensitive to seasonal forcing (R_eff ~ 1). `round()` to whole agents; R_seed
  # takes the remainder. I starts at 0 (sparked by imports). Agents are laid out
  # node-by-node, so node k owns [offset[k], offset[k] + pops[k]); set its first
  # S_seed[k] to S and the rest to R.
  states  <- laser_states()                       # c(S=0, E=1, I=2, R=3, D=-1)
  S_seed  <- round(pops / r0)
  R_seed  <- pops - S_seed
  offsets <- cumsum(c(0L, pops[-n_nodes]))         # 0-based start of each node's block
  state0  <- rep(states[["S"]], capacity)          # S everywhere (incl. reserved slots)
  for (k in seq_len(n_nodes)) {
    nr <- R_seed[k]
    if (nr > 0L) state0[offsets[k] + S_seed[k] + seq_len(nr)] <- states[["R"]]
  }
  people$state$set(state0)
  # timer is already all-zero; no infectious agents to seed a recovery clock for.

  # Seed the census at tick 0 (the first n_nodes entries of each buffer); the remaining
  # nticks-1 columns start zero and are filled by the dynamics / carry-forward.
  zeros <- rep(0L, (nticks - 1L) * n_nodes)
  nodes$S$set(c(S_seed,           zeros))
  nodes$I$set(c(integer(n_nodes), zeros))          # I_0 = 0 in every patch
  nodes$R$set(c(R_seed,           zeros))
  nodes$N$set(c(pops,             zeros))          # N_0 = population

  # The importation schedule as three equal-length integer vectors (the Rust kernel
  # scans them for entries matching the current tick).
  sched_tick  <- as.integer(schedule$tick)
  sched_node  <- as.integer(schedule$node)
  sched_count <- as.integer(schedule$count)

  # ── per-tick dynamics (downstream-first), run nticks-1 times ──────────────────────
  # See endemic_sir.R for the full per-step rationale; the only difference here is that
  # calc_foi reads a time-varying `seasonality` column, so the local rate
  # beta*seasonality[t]*I/N rises and falls through the year.
  run <- function() {
    pb <- if (progress) utils::txtProgressBar(min = 0L, max = nticks - 1L, style = 3) else NULL
    on.exit(if (!is.null(pb)) close(pb), add = TRUE)
    update_every <- max(1L, (nticks - 1L) %/% 100L)
    for (tick in seq_len(nticks - 1L)) {
      t0 <- tick - 1L
      carry_forward_states(list(nodes$S, nodes$I, nodes$R), t0, total = nodes$N)
      sir_step(people$state, people$timer, people$nodeid, people$count,
               nodes$I, nodes$R, nodes$recoveries, t0)
      calc_foi(nodes$I, nodes$N, nodes$beta, nodes$seasonality, network, nodes$foi, t0)
      transmission(people$state, people$timer, people$nodeid, people$count,
                   nodes$foi, nodes$S, nodes$I, nodes$incidence, t0,
                   states[["I"]], inf_duration)
      constant_pop_vitals_sir(people$state, people$timer, people$nodeid, people$count,
                              nodes$death_rate, nodes$S, nodes$I, nodes$R,
                              nodes$births, nodes$deaths, t0)
      people$count <- import_infections(
        people$state, people$timer, people$nodeid, people$count,
        nodes$I, nodes$importations, sched_tick, sched_node, sched_count,
        inf_duration, t0)
      if (!is.null(pb) && (tick %% update_every == 0L || tick == nticks - 1L))
        utils::setTxtProgressBar(pb, tick)
    }
  }
  run()

  # Report. The final census is the last column (state nticks-1) of each buffer;
  # `$values()` returns an nticks x n_nodes matrix, so row `nticks` is it.
  final <- nticks
  cat(sprintf(paste0("run_endemic_sir (seasonal): %d patches, %d -> %d agents (%d imported), %d ticks; ",
                     "final S=%d I=%d R=%d N=%d; incidence=%d, recoveries=%d, births=deaths=%d.\n"),
              n_nodes, n, people$count, total_imports, nticks,
              sum(nodes$S$values()[final, ]), sum(nodes$I$values()[final, ]),
              sum(nodes$R$values()[final, ]), sum(nodes$N$values()[final, ]),
              sum(nodes$incidence$values()), sum(nodes$recoveries$values()),
              sum(nodes$births$values())))
  invisible(list(people = people, nodes = nodes))
}

# ── setup ───────────────────────────────────────────────────────────────────────────
# Two equal patches of 500,000 agents each (same scale as endemic_sir.R).
scenario <- data.frame(
  name       = c("patch_a", "patch_b"),
  population = c(500000L, 500000L)
)
n_nodes <- nrow(scenario)

# Spatial coupling: a small symmetric cross-coupling so the two patches exchange a
# little force of infection (1% each way). `matrix(..., byrow = TRUE)` fills row by row;
# row k is the fraction of patch k's FOI exported to each patch, diagonal 0.
network <- matrix(c(0.00, 0.01,
                    0.01, 0.00), nrow = n_nodes, byrow = TRUE)

# Basic reproduction number and infectious-period distribution (mean shape*scale = 8
# ticks). beta is derived as R0 / mean(inf_duration) inside the runner.
r0           <- 3
inf_duration <- dist_gamma(2, 4)

# A relatively high crude death rate (annual deaths per 1000) — turnover resupplies
# susceptibles to sustain endemic dynamics.
cdr <- 20

# Ten years of daily steps.
nticks <- 3650L

# ── seasonal transmission forcing ─────────────────────────────────────────────────
# A gentle annual sinusoid on the transmission MULTIPLIER: beta is scaled by
# 1 + amplitude*sin(2*pi*day/365), so it runs +/-`seasonal_amplitude` around 1 over a
# 365-day year. We use +/-10%: an exploration sweep (recorded in the package CHANGELOG)
# found this the smallest forcing that yields clear, phase-locked ANNUAL waves while the
# per-year trough never collapses to zero (a larger 15-20% amplitude occasionally faded
# out a patch in the off-season). The dramatic response to such a small wobble is the
# point — seeding at S ~ N/R0 leaves R_eff ~ 1, so the system is critically poised.
# `seq_len(nticks) - 1L` is the 0-based day index; `values_map` (in the runner)
# broadcasts this length-nticks vector across both patches.
seasonal_amplitude <- 0.10
seasonality        <- 1 + seasonal_amplitude * sin(2 * pi * (seq_len(nticks) - 1L) / 365)

# Importation schedule: spark 10 new infectious cases into EACH patch every 30 ticks,
# providing a small endemic floor so an off-season trough can't permanently fade out.
# `seq(0, nticks-2, by=30)` are the (0-based) import ticks within the nticks-1 dynamics
# intervals; `expand.grid` forms every (tick, node) pair, to which we attach a count.
import_ticks   <- seq(0L, nticks - 2L, by = 30L)
schedule       <- expand.grid(tick = import_ticks, node = seq_len(n_nodes) - 1L)
schedule$count <- 10L

# ── run ───────────────────────────────────────────────────────────────────────────
# `system.time(expr)` runs `expr` and returns named timings; `[["elapsed"]]` indexes
# the wall-clock seconds field. `progress = TRUE` draws a text progress bar.
timing <- system.time(
  result <- run_endemic_sir(scenario, network, nticks, inf_duration, r0, cdr, schedule,
                            seasonality = seasonality, progress = TRUE))
cat(sprintf("run_endemic_sir (seasonal) completed in %.3f s\n", timing[["elapsed"]]))

# ── plot the seasonal forcing and the annual epidemic response ──────────────────────
# Resolve the running script's directory so the plot lands in examples/output/ no
# matter the working directory (the same convention as the other examples).
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# `$values()` returns each census buffer as an nticks x n_nodes matrix; `rowSums` totals
# over patches. We plot against time in YEARS (tick / 365).
I_total <- rowSums(result$nodes$I$values())
years   <- (seq_len(nticks) - 1L) / 365

# Render to a PNG (a non-interactive Rscript has no on-screen device). Two stacked
# panels (`par(mfrow = c(rows, cols))`): the gentle forcing on top, the pronounced
# annual epidemic response below, sharing the same x-axis so the phase-locking is
# visible. `abline(v = ...)` drops light gridlines at each year boundary.
plot_path <- file.path(out_dir, "endemic_sir_seasonal.png")
grDevices::png(plot_path, width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))

# Top panel: the seasonal transmission multiplier (the +/-10% driver).
plot(years, seasonality, type = "l", col = "darkorange", lwd = 2,
     xlab = "time (years)", ylab = "beta multiplier",
     main = sprintf("Seasonal transmission forcing (+/-%.0f%%)", 100 * seasonal_amplitude))
graphics::abline(v = 0:10, col = "grey85"); graphics::abline(h = 1, col = "grey60", lty = 3)

# Bottom panel: the infectious prevalence — large annual waves locked to the forcing.
plot(years, I_total, type = "l", col = "firebrick", lwd = 2,
     xlab = "time (years)", ylab = "agents",
     main = "Infectious prevalence (summed over patches): annual epidemic waves")
graphics::abline(v = 0:10, col = "grey85")

graphics::par(op)
grDevices::dev.off()
cat(sprintf("wrote seasonal trajectory plot to %s\n", plot_path))

# Endemic + seasonal check: per-year peak and trough of total prevalence, so the run
# visibly oscillates annually without the trough collapsing to zero. `tapply(x, g, f)`
# applies `f` to `x` split by group `g`; here group = year index.
yr      <- (seq_along(I_total) - 1L) %/% 365L
peaks   <- tapply(I_total, yr, max)
troughs <- tapply(I_total, yr, min)
cat("per-year infectious peak:   ", paste(as.integer(peaks),   collapse = ", "), "\n", sep = "")
cat("per-year infectious trough: ", paste(as.integer(troughs), collapse = ", "), "\n", sep = "")
