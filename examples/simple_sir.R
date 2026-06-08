#!/usr/bin/env Rscript

# Simple spatial SIR over the England & Wales measles patches, run with the high-level
# runner run_model(). run_model owns the disease loop (carry_forward -> step_sir ->
# calc_foi -> transmission); CONSTANT-POPULATION vital dynamics are added through its
# step_exit callback, which calls constant_pop_vitals_sir on the model's people/census each
# tick. Because every death is reborn as a susceptible in the same slot, the population
# stays constant — so no extra agent-array capacity is needed and the closed-population
# run_model is sufficient. (For a GROWING population — births that add agents — pass
# run_model a larger `capacity`; see long_run_squash.R / engwal_measles.R.)
#
# Run from anywhere:  Rscript examples/simple_sir.R

library(razer)

# ── run_sir_model: build + run a spatial SIR with constant-pop vital dynamics ─────────
# `scenario` is a data.frame with one row per patch (name, population, latitude, longitude)
# plus optional integer `I` / `R` seed columns. `network` is the N x N coupling matrix.
# `inf_duration` is the infectious-period Distribution. `r0` sets beta = r0 / mean(D).
# `seasonality` is any values_map-broadcastable transmission modifier. `cdr` is the crude
# death rate (annual deaths per 1,000) driving constant-population turnover; 0 disables it.
run_sir_model <- function(scenario, network, nticks, inf_duration, r0, seasonality = 1, cdr = 0) {
  if (!inherits(inf_duration, "Distribution"))
    stop("`inf_duration` must be a Distribution (e.g. dist_constant(7) or dist_gamma(2, 2))")

  # Vital dynamics live in a step_exit callback. `init` adds the per-node daily death-hazard
  # grid (CDR per 1,000/year -> per person/day) and the birth/death flow Columns the kernel
  # writes; `step_exit` runs constant_pop_vitals_sir AFTER transmission each tick (model$tick
  # is the 0-based interval index), keeping the S/I/R census in sync and recording flows.
  vitals_init <- function(model) {
    nn <- model$nodes$count
    model$nodes$death_rate <- values_map(cdr / 1000 / 365, nticks, nn)
    model$nodes$births <- allocate_vector("i32", nticks - 1L, nn)
    model$nodes$deaths <- allocate_vector("i32", nticks - 1L, nn)
  }
  vitals_step <- function(model) {
    constant_pop_vitals_sir(model$people$state, model$people$timer, model$people$nodeid,
      model$people$count, model$nodes$death_rate, model$nodes$S, model$nodes$I, model$nodes$R,
      model$nodes$births, model$nodes$deaths, model$tick)
  }

  # run_model builds the agent population from the scenario (seeding I with an infectious
  # timer and R as immune), advances the SIR kernels in the correct order, and returns the
  # `model` environment. seed = 1L makes the stochastic run reproducible.
  m <- run_model(scenario = scenario, model = "SIR", nticks = nticks, r0 = r0,
                 infectious_period = inf_duration, network = network,
                 seasonality = seasonality, seed = 1L,
                 init      = if (cdr > 0) vitals_init else NULL,
                 step_exit = if (cdr > 0) vitals_step else NULL)

  # Report. The final census is the last recorded tick of each S/I/R Column.
  nn <- m$nodes$count
  fin <- function(col) sum(col$values()[nticks, ])
  born <- if (!is.null(m$nodes$births)) sum(m$nodes$births$values()) else 0L
  cat(sprintf(paste0("run_sir_model: %d patches, %d agents, %d ticks; seeded I=%d R=%d; ",
                     "final S=%d I=%d R=%d; incidence=%d, recoveries=%d, births=deaths=%d.\n"),
              nn, m$people$count, nticks,
              if ("I" %in% names(scenario)) sum(scenario$I) else 0L,
              if ("R" %in% names(scenario)) sum(scenario$R) else 0L,
              fin(m$nodes$S), fin(m$nodes$I), fin(m$nodes$R),
              sum(m$nodes$incidence$values()), sum(m$nodes$recovery$values()), born))
  invisible(m)
}

# ── setup ───────────────────────────────────────────────────────────────────────
# Locate this script's directory so the data file resolves no matter the working directory.
args       <- commandArgs(trailingOnly = FALSE)        # all launch args incl. --file=
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."

# England & Wales measles patches (954 registration districts): one row per node with its
# name, initial (1944) population, and latitude/longitude.
scenario <- read.csv(file.path(script_dir, "data", "EnglandWalesMeasles_places.csv"))

# Pairwise great-circle distances (km) between patches, then a radiation-model coupling
# network (Simini et al., Nature 2012). `row_normalizer` caps each patch's total emigration
# fraction (here 10%) so the network rows are valid export fractions for calc_foi.
distance_matrix <- distances(scenario$latitude, scenario$longitude)
network <- radiation(scenario$population, distance_matrix, k = 1, include_home = FALSE)
network <- row_normalizer(network, max_rowsum = 0.1)

# Initial conditions as integer I / R columns on the scenario: 5 infectious per patch and
# ~5% immune, capped so I + R never exceeds a patch's population.
scenario$I <- pmin(5L, scenario$population)
scenario$R <- pmin(as.integer(0.05 * scenario$population), scenario$population - scenario$I)

inf_duration <- dist_gamma(2, 2)        # infectious period (mean shape*scale = 4 ticks)
R0           <- 2                       # beta = R0 / mean(inf_duration)
nticks       <- 10L                     # daily time steps
# Seasonal transmission modifier: a per-tick sinusoid (values_map also accepts a scalar,
# a per-node vector, or a full nticks x n_nodes matrix).
seasonality  <- 1 + 0.3 * sin(2 * pi * (seq_len(nticks) - 1L) / nticks)
cdr          <- 20                      # crude death rate (annual per 1,000) for vitals

# ── run ───────────────────────────────────────────────────────────────────────
timing <- system.time(run_sir_model(scenario, network, nticks, inf_duration, R0, seasonality, cdr))
cat(sprintf("run_sir_model completed in %.3f s\n", timing[["elapsed"]]))
