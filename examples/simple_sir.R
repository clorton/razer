#!/usr/bin/env Rscript

# Simple spatial SIR example over the England & Wales measles patches. This script
# defines a self-contained run_sir_model() that wires the Column per-tick kernels
# (carry_forward_states -> calc_foi -> step_sir -> transmission -> vital dynamics),
# does some setup, then calls it.
#
# Run from anywhere:  Rscript examples/simple_sir.R

# `library(pkg)` attaches a package so its exported names resolve unqualified
# (like a wildcard `import`); the package name is given bare, not as a string.
library(razer)

# ── run_sir_model: assemble parameters and run the Column SIR loop ──────────────
# `function(args) body` is a first-class closure value, bound here with `<-`.
# `scenario` is a data.frame with one row per patch (name, population, latitude,
# longitude); it is the model's geographic scaffold. Optional integer `I` / `R`
# columns give the number of agents to start infectious / recovered in each patch.
# `network` is the N x N spatial coupling matrix (fraction of each patch's force of
# infection exported to every other patch). `nticks` is the number of daily time
# steps to simulate. `inf_duration` is a Distribution (e.g. dist_gamma(2, 2)) from
# which each infectious agent's recovery timer is drawn. `r0` is the basic
# reproduction number; the transmission coefficient is `beta = r0 / mean(inf_duration)`.
# `seasonality` is a transmission modifier accepted in any broadcastable form (a
# scalar, a per-node or per-tick vector, or a full nticks x nnodes matrix). `cdr` is
# the crude death rate (annual deaths per 1000 people) driving constant-population
# vital dynamics; 0 (the default) disables them.
run_sir_model <- function(scenario, network, nticks, inf_duration, r0, seasonality = 1, cdr = 0) {
  # `inherits(x, "Distribution")` checks the S3 class of the extendr handle.
  if (!inherits(inf_duration, "Distribution"))
    stop("`inf_duration` must be a Distribution (e.g. dist_constant(7) or dist_gamma(2, 2))")

  # beta from R0: R0 = beta * mean infectious duration, so beta = R0 / mean. Estimate
  # the mean by sampling the distribution (a large batch; `mean()` averages it).
  beta <- r0 / mean(inf_duration$sample_n(100000L))

  # Total agent count = sum of every patch's population. `sum()` over an integer
  # column returns an integer scalar. `n_nodes` is the patch count.
  n       <- sum(scenario$population)
  n_nodes <- nrow(scenario)

  # `new.env()` makes a fresh ENVIRONMENT: unlike ordinary R values (which are
  # copy-on-modify), an environment has reference semantics — assigning into it
  # mutates the same object everyone holds, like a mutable struct / object. We use
  # it as a lightweight "people" record whose named entries are agent-property
  # arrays. `env$name <- value` (the `$<-` method) assigns a member by name.
  people <- new.env()
  people$count    <- n   # number of live agents
  people$capacity <- n   # allocated array length (room to grow); equals count for now
  # `allocate_scalar()` (Rust) returns an opaque Column handle backed by a
  # Rust-owned, zero-filled uint8 buffer of length n, holding each agent's disease
  # state. The step kernels mutate this buffer in place (no copies); `$dtype()`
  # and `$length()` query it, `$values()` copies a snapshot back to R.
  people$state    <- allocate_scalar("u8", n)

  # Per-agent node id — which patch each agent lives in — as a compact uint16.
  # nodeid is 0-BASED (0..N-1), matching the Rust step kernels and the rest of
  # razer: the kernels index per-node arrays (infectious counts, populations, the
  # coupling network) DIRECTLY by this value with no offset, so storing it 1-based
  # would cost a `- 1` on every hot-path access. R-side joins back to `scenario`
  # rows add 1 instead (e.g. scenario[nodeid + 1, ]) — a rare, cheap conversion.
  people$nodeid   <- allocate_scalar("u16", n)
  # `rep(x, times = counts)` repeats each node id by that patch's population, so
  # the first population[1] agents are node 0, the next population[2] are node 1,
  # and so on. `seq_len(N) - 1L` is the 0-based id vector. `$set()` copies it once
  # into the Rust buffer.
  people$nodeid$set(rep(seq_len(nrow(scenario)) - 1L, times = scenario$population))

  # Per-agent countdown timer (uint16), zero-initialized. The step kernels set it
  # when an agent enters a timed compartment (e.g. the infectious period on S->I)
  # and decrement it each tick; a uint16 holds durations up to 65535 ticks (long
  # immunity / maternal periods need more than a uint8's 255).
  people$timer    <- allocate_scalar("u16", n)

  # `nodes` is the per-patch record. `nticks` is the number of recorded time points
  # (states 0..nticks-1): state 0 is the seeded initial condition and state nticks-1
  # is the final state. The dynamics therefore run nticks-1 times (one per interval
  # between consecutive states); running them on the last state would write a state
  # nticks that falls outside the simulation window. So the S/I/R CENSUS buffers are
  # nticks x n_nodes (maintained incrementally: each step carries column t forward to
  # t+1 and applies only the deltas), and the `foi` / `recoveries` / `incidence`
  # per-interval FLOW reports are (nticks-1) x n_nodes. All are time-major (each
  # tick's per-node row is contiguous).
  nodes <- new.env()
  nodes$count      <- n_nodes
  nodes$S          <- allocate_vector("i32", nticks,      n_nodes)
  nodes$I          <- allocate_vector("i32", nticks,      n_nodes)
  nodes$R          <- allocate_vector("i32", nticks,      n_nodes)
  nodes$foi        <- allocate_vector("f64", nticks - 1L, n_nodes)
  nodes$recoveries <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$incidence  <- allocate_vector("i32", nticks - 1L, n_nodes)
  # Per-node population (the FOI denominator): constant here, but kept as a full
  # time x node grid because vital dynamics will eventually make it vary by tick.
  nodes$N          <- values_map(scenario$population, nticks, n_nodes)
  # Transmission inputs as time x node grids (built by values_map from whatever
  # shape was supplied): the beta coefficient and the seasonal modifier. calc_foi
  # multiplies them into the local force of infection.
  nodes$beta        <- values_map(beta,        nticks, n_nodes)
  nodes$seasonality <- values_map(seasonality, nticks, n_nodes)
  # Vital dynamics: per-node daily death-HAZARD rate from the crude death rate
  # (annual deaths per 1000 -> daily per person), plus per-interval birth/death flows.
  nodes$death_rate  <- values_map(cdr / 1000 / 365, nticks, n_nodes)
  nodes$births      <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$deaths      <- allocate_vector("i32", nticks - 1L, n_nodes)

  # ── seed initial infections / immunity from the scenario ──────────────────────
  # Agents are laid out node-by-node (see nodeid), so node k owns the contiguous
  # block [offset[k], offset[k] + population[k]). For each node we move its first
  # I[k] agents to the I state and the next R[k] to R. `%in%` tests for a column;
  # a missing I/R column means zero seeding for that state.
  states  <- laser_states()                          # c(S=0, E=1, I=2, R=3, D=-1)
  pops    <- scenario$population
  offsets <- cumsum(c(0L, pops[-n_nodes]))           # 0-based start of each node's block
  I_seed  <- if ("I" %in% names(scenario)) as.integer(scenario$I) else integer(n_nodes)
  R_seed  <- if ("R" %in% names(scenario)) as.integer(scenario$R) else integer(n_nodes)
  if (any(I_seed < 0L) || any(R_seed < 0L))
    stop("`scenario$I` and `scenario$R` must be non-negative")
  # I + R must fit within each patch's population, else there aren't enough agents
  # to place into both states — an irreconcilable initial condition. `which()`
  # returns the offending row indices; name the first for a actionable message.
  bad <- which(I_seed + R_seed > pops)
  if (length(bad) > 0L)
    stop(sprintf(
      "seeded I + R exceeds population in %d patch(es); e.g. patch %d: I=%d + R=%d > population=%d",
      length(bad), bad[1L], I_seed[bad[1L]], R_seed[bad[1L]], pops[bad[1L]]))

  # Build the initial state vector (S everywhere, then per-node I then R) and write
  # it into the Rust buffer with one `$set()`.
  state0 <- rep(states[["S"]], n)
  for (k in seq_len(n_nodes)) {
    base <- offsets[k]; ni <- I_seed[k]; nr <- R_seed[k]
    if (ni > 0L) state0[base + seq_len(ni)]      <- states[["I"]]
    if (nr > 0L) state0[base + ni + seq_len(nr)] <- states[["R"]]
  }
  people$state$set(state0)

  # Each infectious agent gets a recovery timer drawn from `inf_duration`, rounded
  # to whole ticks and clamped to a u8's [1, 255]; everyone else stays at 0.
  timer0  <- rep(0L, n)
  total_I <- sum(I_seed)
  if (total_I > 0L) {
    draws <- pmax(1L, pmin(255L, as.integer(round(inf_duration$sample_n(total_I)))))
    pos <- 0L
    for (k in seq_len(n_nodes)) {
      ni <- I_seed[k]
      if (ni > 0L) {
        timer0[offsets[k] + seq_len(ni)] <- draws[pos + seq_len(ni)]
        pos <- pos + ni
      }
    }
  }
  people$timer$set(timer0)

  # Seed the census at tick 0 (S = population - I - R per node). `$set()` writes the
  # whole buffer; column 0 is the first n_nodes entries, the remaining nticks-1
  # columns start at 0.
  nodes$S$set(c(pops - I_seed - R_seed, rep(0L, (nticks - 1L) * n_nodes)))
  nodes$I$set(c(I_seed,                 rep(0L, (nticks - 1L) * n_nodes)))
  nodes$R$set(c(R_seed,                 rep(0L, (nticks - 1L) * n_nodes)))

  # `run` advances the simulation: each tick runs the per-tick step kernels in
  # downstream-first order. `tick - 1L` converts R's 1-based loop counter to the
  # 0-based tick index the kernels expect.
  #   carry_forward:    seed tick t+1 of each census counter with tick t's value
  #                     (an SEIR model would also carry E; users can carry e.g. V).
  #   calc_foi:         FOI[t] from the infectious census I[t+1] BEFORE recovery, so an
  #                     agent counts on its recovery tick too. With direct S->I (no E),
  #                     the agent enters I *after* this tally, so it is not counted on
  #                     its entry tick; counting it on its recovery tick instead gives
  #                     the full infectious period D — i.e. R0 = beta * D, not D - 1.
  #                     (Run calc_foi BEFORE step_sir; this is the SIR ordering.)
  #   step_sir:         recover (I -> R), updating the carried-forward census at t+1.
  #   transmission:     infect susceptibles from FOI[t] (S -> I) at t+1. For SIR the
  #                     receiving state is I directly (an SEIR model would pass E and
  #                     an incubation-period distribution instead).
  #   vital dynamics:   constant-population births/deaths (deaths reborn susceptible),
  #                     keeping the S/I/R census at t+1 in sync and recording flows.
  # Run dynamics nticks-1 times (the intervals between the nticks recorded states);
  # the final state nticks-1 has nothing transitioning out of the window.
  run <- function() {
    for (tick in seq_len(nticks - 1L)) {
      t0 <- tick - 1L
      carry_forward(nodes$S, t0)
      carry_forward(nodes$I, t0)
      carry_forward(nodes$R, t0)
      calc_foi(nodes$I, nodes$N, nodes$beta, nodes$seasonality, network, nodes$foi, t0)
      # step_sir(absorbing = R): I->R recovery. SIR has no M/E, so waned/onset are 0.
      rec <- step_sir(people$state, people$timer, people$nodeid, people$count,
                      nodes$count, inf_duration, states[["R"]])
      move_count(nodes$I, nodes$R, rec$cleared, t0)
      nodes$recoveries$set_col(t0, rec$cleared)
      # transmission S->I returns new infections per node; apply the S->I delta.
      inf <- transmission(people$state, people$timer, people$nodeid, people$count,
                          nodes$foi, t0, states[["I"]], inf_duration)
      move_count(nodes$S, nodes$I, inf, t0)
      nodes$incidence$set_col(t0, inf)
      constant_pop_vitals_sir(people$state, people$timer, people$nodeid, people$count,
                              nodes$death_rate, nodes$S, nodes$I, nodes$R,
                              nodes$births, nodes$deaths, t0)
    }
  }
  run()

  # Report success. The final census is the last column (state nticks-1) of S/I/R;
  # `$values()` returns an nticks x n_nodes matrix, so row `nticks` is it.
  final <- nticks
  cat(sprintf(paste0("run_sir_model: %d patches, %d agents, %d ticks; seeded I=%d R=%d; ",
                     "final S=%d I=%d R=%d; incidence=%d, recoveries=%d, births=deaths=%d.\n"),
              n_nodes, people$count, nticks, sum(I_seed), sum(R_seed),
              sum(nodes$S$values()[final, ]), sum(nodes$I$values()[final, ]),
              sum(nodes$R$values()[final, ]),
              sum(nodes$incidence$values()), sum(nodes$recoveries$values()),
              sum(nodes$births$values())))
  invisible(list(people = people, nodes = nodes))
}

# ── setup ───────────────────────────────────────────────────────────────────────
# Locate this script's directory so the data file resolves no matter the working
# directory. `commandArgs(FALSE)` returns ALL launch args including Rscript's own
# `--file=<path>`; `grep(..., value = TRUE)` returns the matching element, `sub`
# strips the prefix, and `dirname` takes the directory part. The `else "."`
# fallback covers being sourced interactively.
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."

# Load the England & Wales measles patches (954 registration districts): one row
# per node with its name, initial (1944) population, and latitude/longitude.
# `read.csv` returns a data.frame; `file.path` joins paths OS-correctly. The CSV
# is produced from EnglandWalesMeasles.py by examples/data/convert_measles.py.
scenario <- read.csv(file.path(script_dir, "data", "EnglandWalesMeasles_places.csv"))

# Pairwise great-circle distances (km) between every patch, from their
# latitude/longitude. `distances()` returns a symmetric N x N matrix with a zero
# diagonal; it is the geographic input from which the spatial coupling network is
# built (e.g. via a gravity model). `$` reaches a data.frame column by name.
distance_matrix <- distances(scenario$latitude, scenario$longitude)

# Spatial coupling network from the radiation model (Simini et al., Nature 2012):
# given each patch's population and the pairwise distances, it estimates the
# migration weight between every ordered pair of patches. `include_home = FALSE`
# excludes a patch's own population from the "intervening opportunities" sum.
network <- radiation(scenario$population, distance_matrix, k = 1, include_home = FALSE)

# The transmission kernel reads the network as the FRACTION of each patch's force
# of infection exported to every other patch, so a row sum is that patch's total
# emigration fraction and must be at most 1. `row_normalizer` proportionally
# scales down any row that exceeds the cap; here we limit every patch to at most
# 10% total emigration.
network <- row_normalizer(network, max_rowsum = 0.1)

# Initial conditions, attached to `scenario` as integer `I` / `R` columns: a few
# infectious agents per patch and ~5% initial immunity. `pmin` caps the counts so
# `I + R` never exceeds a patch's population. `as.integer` truncates toward zero.
scenario$I <- pmin(5L, scenario$population)                       # 5 infectious / patch
scenario$R <- pmin(as.integer(0.05 * scenario$population),        # ~5% recovered/immune
                   scenario$population - scenario$I)

# Infectious-period distribution: each infected agent's recovery timer is drawn
# from this (mean shape*scale = 4 ticks). Drives the recovery clock.
inf_duration <- dist_gamma(2, 2)

# Basic reproduction number; run_sir_model derives beta = R0 / mean(inf_duration).
R0 <- 2

# Number of daily time steps to simulate.
nticks <- 10L

# Seasonal transmission modifier: here a per-tick value (length nticks) — a gentle
# sinusoid for illustration. values_map() would equally accept a scalar, a
# per-node vector, or a full nticks x nnodes matrix.
seasonality <- 1 + 0.3 * sin(2 * pi * (seq_len(nticks) - 1L) / nticks)

# Crude death rate (annual deaths per 1000 people) for constant-population vital
# dynamics; births equal deaths, resupplying susceptibles for endemic dynamics.
cdr <- 20

# ── run ───────────────────────────────────────────────────────────────────────
# `system.time(expr)` runs `expr` and returns a named numeric vector of timings;
# `[["elapsed"]]` name-indexes the wall-clock seconds field.
timing <- system.time(run_sir_model(scenario, network, nticks, inf_duration, R0, seasonality, cdr))
cat(sprintf("run_sir_model completed in %.3f s\n", timing[["elapsed"]]))
