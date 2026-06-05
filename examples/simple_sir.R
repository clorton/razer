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
# to every other patch).
run_sir_model <- function(scenario, network) {
  # Total agent count = sum of every patch's population. `sum()` over an integer
  # column returns an integer scalar.
  n <- sum(scenario$population)

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

  # TODO: assign agents to patches, seed infections, build SIR parameters, run.

  # Report success. `cat()` writes to stdout; `sprintf` is C-style formatting.
  # `invisible()` returns its argument WITHOUT auto-printing at the top level.
  cat(sprintf("run_sir_model: %d patches, network %d x %d; people$state = %s[%d] (capacity %d).\n",
              nrow(scenario), nrow(network), ncol(network),
              people$state$dtype(), people$state$length(), people$capacity))
  invisible(people)
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

# ── run ───────────────────────────────────────────────────────────────────────
# `system.time(expr)` runs `expr` and returns a named numeric vector of timings;
# `[["elapsed"]]` name-indexes the wall-clock seconds field.
timing <- system.time(run_sir_model(scenario, network))
cat(sprintf("run_sir_model completed in %.3f s\n", timing[["elapsed"]]))
