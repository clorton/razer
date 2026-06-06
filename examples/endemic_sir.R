#!/usr/bin/env Rscript

# Endemic SIR example: a small two-patch metapopulation driven to ENDEMICITY by
# constant-population vital dynamics (a high crude death rate continually resupplies
# susceptibles) plus periodic IMPORTATIONS that re-spark transmission so a stochastic
# fade-out doesn't end the epidemic. It exercises three pieces beyond simple_sir.R:
#   * import_infections()      — activates RESERVED agent slots as new infectious cases
#                                from a schedule (so `people` is allocated with capacity
#                                ABOVE the initial count to leave room for imports);
#   * carry_forward_states()   — carries the S/I/R census forward AND totals it into N
#                                (the population / FOI denominator) in one call;
#   * constant_pop_vitals_sir()— births = deaths, reborn susceptible, keeping S+I+R=N.
#
# Run from anywhere:  Rscript examples/endemic_sir.R

# `library(pkg)` attaches a package so its exported names resolve unqualified (like a
# wildcard `import`); the package name is given bare, not as a string.
library(razer)

# ── run_endemic_sir: build the model, run it, report ────────────────────────────────
# `function(args) body` is a first-class closure value bound with `<-`. `scenario` is a
# data.frame with one row per patch (name, population); `network` is the N x N spatial
# coupling matrix; `nticks` the number of recorded daily states; `inf_duration` the
# infectious-period Distribution; `r0` the basic reproduction number (beta = r0/mean);
# `cdr` the crude death rate (annual deaths per 1000); `schedule` a data.frame of
# importations with integer columns tick (0-based), node (0-based), count.
run_endemic_sir <- function(scenario, network, nticks, inf_duration, r0, cdr, schedule) {
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
  # The S/I/R CENSUS buffers are nticks x n_nodes (state 0 is the initial condition,
  # state nticks-1 the final; maintained incrementally — each tick carries column t
  # forward to t+1 and applies only deltas). N (population / FOI denominator) is the
  # same shape, recomputed each tick as S+I+R so it tracks the agents added by imports.
  # The per-interval FLOW reports are (nticks-1) x n_nodes. All are time-major (each
  # tick's per-node row is contiguous).
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
  # values_map from a scalar — here no spatial or temporal variation. The death rate
  # is the per-node daily death HAZARD (crude death rate per 1000 per year -> daily).
  nodes$beta         <- values_map(beta,             nticks, n_nodes)
  nodes$seasonality  <- values_map(1,                nticks, n_nodes)
  nodes$death_rate   <- values_map(cdr / 1000 / 365, nticks, n_nodes)

  # ── seed the initial condition: 1/R0 susceptible, the rest recovered, I = 0 ───────
  # At the endemic equilibrium of an SIR model the susceptible fraction is 1/R0, so we
  # start there (the rest immune) and let the imports + vital turnover sustain it.
  # `round()` to whole agents; R_seed takes the remainder. I starts at 0 (sparked by
  # imports). Agents are laid out node-by-node, so node k owns the contiguous block
  # [offset[k], offset[k] + pops[k]); set its first S_seed[k] to S and the rest to R.
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
  # `tick - 1L` converts R's 1-based loop counter to the 0-based tick index the kernels
  # expect (t0). Order:
  #   carry_forward_states: copy each census column t0 -> t0+1 and set N[t0+1] = S+I+R
  #                         (one call does both the carry and the population total).
  #   sir_step:             recover expired infectious (I -> R) at column t0+1.
  #   calc_foi:             FOI[t0] from the post-recovery census I[t0+1] / N[t0+1].
  #   transmission:         infect susceptibles (S -> I) at t0+1, timer from inf_duration.
  #   constant_pop_vitals_sir: deaths (reborn susceptible) = births, census kept in sync.
  #   import_infections:    activate reserved slots as new infectious cases per the
  #                         schedule; returns the grown live count. Placed last so the
  #                         imports seed the NEXT tick's force of infection (the same
  #                         one-tick delay new S->I infections already have), and so the
  #                         next carry_forward_states folds them into N.
  run <- function() {
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
    }
  }
  run()

  # Report. The final census is the last column (state nticks-1) of each buffer;
  # `$values()` returns an nticks x n_nodes matrix, so row `nticks` is it.
  final <- nticks
  cat(sprintf(paste0("run_endemic_sir: %d patches, %d -> %d agents (%d imported), %d ticks; ",
                     "final S=%d I=%d R=%d N=%d; incidence=%d, recoveries=%d, births=deaths=%d.\n"),
              n_nodes, n, people$count, total_imports, nticks,
              sum(nodes$S$values()[final, ]), sum(nodes$I$values()[final, ]),
              sum(nodes$R$values()[final, ]), sum(nodes$N$values()[final, ]),
              sum(nodes$incidence$values()), sum(nodes$recoveries$values()),
              sum(nodes$births$values())))
  invisible(list(people = people, nodes = nodes))
}

# ── setup ───────────────────────────────────────────────────────────────────────────
# Two equal patches of 500,000 agents each. `data.frame` builds the scenario table.
scenario <- data.frame(
  name       = c("patch_a", "patch_b"),
  population = c(500000L, 500000L)
)
n_nodes <- nrow(scenario)

# Spatial coupling: a small symmetric cross-coupling so the two patches exchange a
# little force of infection (1% each way). `matrix(..., byrow = TRUE)` fills row by row;
# row k is the fraction of patch k's FOI exported to each patch, diagonal 0 (a patch's
# own contribution is implicit).
network <- matrix(c(0.00, 0.01,
                    0.01, 0.00), nrow = n_nodes, byrow = TRUE)

# Basic reproduction number and infectious-period distribution (mean shape*scale = 8
# ticks). beta is derived as R0 / mean(inf_duration) inside the runner.
r0           <- 3
inf_duration <- dist_gamma(2, 4)

# A relatively high crude death rate (annual deaths per 1000) with NO spatial or
# temporal variation — turnover resupplies susceptibles to sustain endemic dynamics.
cdr <- 20

# Two years of daily steps.
nticks <- 730L

# Importation schedule: spark 10 new infectious cases into EACH patch every 30 ticks,
# keeping the epidemic from stochastically fading out. `seq(0, nticks-2, by=30)` are the
# (0-based) import ticks within the nticks-1 dynamics intervals; `expand.grid` forms
# every (tick, node) pair, to which we attach a constant count.
import_ticks   <- seq(0L, nticks - 2L, by = 30L)
schedule       <- expand.grid(tick = import_ticks, node = seq_len(n_nodes) - 1L)
schedule$count <- 10L

# ── run ───────────────────────────────────────────────────────────────────────────
# `system.time(expr)` runs `expr` and returns named timings; `[["elapsed"]]` indexes
# the wall-clock seconds field.
timing <- system.time(
  result <- run_endemic_sir(scenario, network, nticks, inf_duration, r0, cdr, schedule))
cat(sprintf("run_endemic_sir completed in %.3f s\n", timing[["elapsed"]]))

# Endemic check: print the I prevalence trajectory (summed over patches) at a few
# checkpoints so the run visibly settles around a non-zero level rather than going
# extinct or exploding. `seq(...)` picks evenly-spaced recorded states; `rowSums` totals
# the per-patch matrix over its columns.
I_total     <- rowSums(result$nodes$I$values())
checkpoints <- unique(round(seq(1, nticks, length.out = 6)))
cat("infectious prevalence (total over patches) at ticks ",
    paste(checkpoints - 1L, collapse = ", "), ":\n  ",
    paste(I_total[checkpoints], collapse = ", "), "\n", sep = "")
