# Tests for run_model(): the high-level runner for the closed-population menagerie. It
# builds people/nodes, seeds, runs the correct per-tick loop, and applies the census
# deltas. Closed population (no births/deaths), so the living compartments sum to N at
# every tick. Written given-when-then.

# Run each model with the durations it needs; returns the result.
run_one <- function(model, nticks = 60L, n = 20000L, seed = 1L) {
  scenario <- data.frame(population = n, I = 50L)
  run_model(scenario, model, nticks = nticks, r0 = 2.5,
            infectious_period = 8, incubation_period = 5, immunity_period = 60, seed = seed)
}

# Sum the present census compartments at each recorded tick (an n_ticks-length vector).
living_by_tick <- function(nodes, model) {
  cols <- list(nodes$S, nodes$I)
  if (grepl("E", model)) cols <- c(cols, list(nodes$E))
  if (grepl("R", model)) cols <- c(cols, list(nodes$R))
  Reduce(`+`, lapply(cols, function(c) rowSums(c$values())))
}

test_that("run_model runs every model and conserves the closed population", {
  # Given each of the eight models
  # When run for 60 ticks at R0 = 2.5
  # Then the present compartments sum to N at every recorded tick (no births/deaths), and
  #      some transmission has occurred (final S < initial S). Failure would mean a
  #      mis-applied census delta (the desync the runner exists to prevent) or a model
  #      mis-wired.
  for (model in c("SI", "SEI", "SIS", "SEIS", "SIR", "SEIR", "SIRS", "SEIRS")) {
    res <- run_one(model)
    living <- living_by_tick(res$nodes, model)
    expect_true(all(living == 20000L), info = paste(model, "conserves population"))
    s <- rowSums(res$nodes$S$values())
    expect_lt(s[length(s)], s[1L])             # susceptibles fell -> transmission happened
  }
})

test_that("run_model is reproducible under a seed", {
  # Given the same seed
  # When an SEIR model is run twice
  # Then the infectious trajectories are identical. Failure would mean the runner does not
  #      thread the seed through deterministically.
  a <- run_one("SEIR", seed = 7L)
  b <- run_one("SEIR", seed = 7L)
  expect_identical(a$nodes$I$values(), b$nodes$I$values())
})

test_that("SIS returns to susceptible (no R compartment); SIR accumulates R", {
  # SIS: I clears back to S, so there is no R census and S recovers after the peak.
  sis <- run_one("SIS")
  expect_null(sis$nodes$R)
  # SIR: recovereds accumulate monotonically (lifelong immunity, closed population).
  sir <- run_one("SIR")
  r <- rowSums(sir$nodes$R$values())
  expect_true(all(diff(r) >= 0))
  expect_gt(r[length(r)], 0)
})

test_that("run_model validates its inputs", {
  scen <- data.frame(population = 1000L, I = 10L)
  expect_error(run_model(scen, "SXYZ", nticks = 10L, r0 = 2, infectious_period = 5), "must be one of")
  expect_error(run_model(scen, "SEIR", nticks = 10L, r0 = 2, infectious_period = 5), "incubation_period")
  expect_error(run_model(scen, "SIRS", nticks = 10L, r0 = 2, infectious_period = 5), "immunity_period")
  expect_error(run_model(scen, "SIR", nticks = 1L, r0 = 2, infectious_period = 5), "nticks")
})
