#!/usr/bin/env Rscript

# Simple SIR example built on the high-level run_sir() runner (R/models.R),
# rather than wiring the per-tick step kernels together by hand. This script
# defines run_sir_model(), does some setup, then calls it.
#
# Run from anywhere:  Rscript examples/simple_sir.R

# `library(pkg)` attaches a package so its exported names resolve unqualified
# (like a wildcard `import`); the package name is given bare, not as a string.
library(razer)

# ── run_sir_model: assemble parameters and hand off to run_sir() ────────────────
# `function(args) body` is a first-class closure value, bound here with `<-`.
# `scenario` is a data.frame with one row per patch (name, population, latitude,
# longitude); it is the model's geographic scaffold. `network` is the N x N
# spatial coupling matrix (fraction of each patch's force of infection exported
# to every other patch). `nticks` is the number of daily time steps to simulate.
run_sir_model <- function(scenario, network, nticks) {
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

  # Per-agent countdown timer (uint8), zero-initialized. The step kernels set it
  # when an agent enters a timed compartment (e.g. the infectious period on S->I)
  # and decrement it each tick; a uint8 holds durations up to 255 ticks.
  people$timer    <- allocate_scalar("u8", n)

  # `nodes` is the per-patch record holding per-tick REPORTS. `recoveries` is a
  # 2-D buffer of shape (nticks, n_nodes) — time first, node second — laid out so
  # each tick's per-node row is contiguous in memory; `sir_step` tallies recoveries
  # into the current tick's slice.
  nodes <- new.env()
  nodes$count      <- n_nodes
  nodes$recoveries <- allocate_vector("u32", nticks, n_nodes)

  # TODO: seed infections, build the rest of the SIR parameters.

  # `run` advances the simulation: each tick runs the per-tick step kernels (just
  # the SIR recovery step for now; transmission etc. will join it here). `tick - 1L`
  # converts R's 1-based loop counter to the 0-based time-slice index sir_step expects.
  run <- function() {
    for (tick in seq_len(nticks)) {
      sir_step(people$state, people$timer, people$nodeid, people$count,
               nodes$recoveries, tick - 1L)
    }
  }
  run()

  # Report success. `cat()` writes to stdout; `sprintf` is C-style formatting.
  # `sum(nodes$recoveries$values())` totals the recovery report. `invisible()`
  # returns its argument WITHOUT auto-printing at the top level.
  cat(sprintf("run_sir_model: %d patches, %d ticks; people = {state %s, nodeid %s, timer %s}[%d]; total recoveries = %d.\n",
              n_nodes, nticks,
              people$state$dtype(), people$nodeid$dtype(), people$timer$dtype(),
              people$count, sum(nodes$recoveries$values())))
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

# Number of daily time steps to simulate (placeholder until the full SIR loop and
# seeding are wired up).
nticks <- 10L

# ── run ───────────────────────────────────────────────────────────────────────
# `system.time(expr)` runs `expr` and returns a named numeric vector of timings;
# `[["elapsed"]]` name-indexes the wall-clock seconds field.
timing <- system.time(run_sir_model(scenario, network, nticks))
cat(sprintf("run_sir_model completed in %.3f s\n", timing[["elapsed"]]))
