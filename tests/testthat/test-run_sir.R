# Tests for run_sir(): the agent-based SIR runner that builds a population from a
# node data.frame, runs the downstream-first SIR loop, and returns the people
# LaserFrame with per-node trajectories (S/I/R) and flows (incidence/recovery)
# attached via attr(model, "report").
#
# Written given-when-then. testthat idioms: `expect_*` assertions; `L` marks an
# integer literal; `set.seed` fixes the RNG so stochastic runs are reproducible.

# Helper: a deterministic two-node scenario reused across tests.
two_nodes <- function() data.frame(population = c(10000L, 5000L))

test_that("run_sir assigns agents to nodes by population count", {
  # Given a two-node table with populations 10000 and 5000
  # When the model is built and run
  # Then the people frame holds 15000 agents and the per-node membership counts
  #      match the requested populations exactly.
  # Failure would mean the rep(node, times = population) assignment is wrong and
  # every downstream per-node tally would be misattributed.
  nodes <- two_nodes()
  set.seed(1L)
  model <- run_sir(nodes, beta = 0.3, infectious_period = 7, nticks = 5L,
                   seed = c(5L, 0L))

  expect_equal(model$count, 15000L)
  node_counts <- tabulate(model$node + 1L, nbins = 2L)  # +1L: 0-based -> 1-based
  expect_equal(node_counts, c(10000L, 5000L))
})

test_that("run_sir seeds the requested number of infectious agents per node", {
  # Given a per-node seed of c(10, 0)
  # When the model is built
  # Then the initial (column 1) infectious counts are exactly 10 and 0, and the
  #      single susceptible-only initial state holds elsewhere.
  # Failure would indicate the seeding block placed infections in the wrong
  #      agents or the wrong node block.
  nodes <- two_nodes()
  set.seed(2L)
  model <- run_sir(nodes, beta = 0.3, infectious_period = 7, nticks = 10L,
                   seed = c(10L, 0L))
  report <- attr(model, "report")

  expect_equal(report$I[, 1L], c(10L, 0L))
  expect_equal(report$S[, 1L], c(9990L, 5000L))
  expect_equal(report$R[, 1L], c(0L, 0L))
})

test_that("run_sir conserves population within each node at every tick", {
  # Given any SIR run (no births or deaths)
  # When trajectories are recorded over the whole horizon
  # Then S + I + R equals the node population in every node at every tick.
  # Failure would mean an agent was lost or double-counted by the tally or a
  #      kernel — a correctness red flag for the whole pipeline.
  nodes <- two_nodes()
  set.seed(3L)
  model <- run_sir(nodes, beta = 0.4, infectious_period = 6, nticks = 80L,
                   seed = c(20L, 5L))
  report <- attr(model, "report")

  totals <- report$S + report$I + report$R   # element-wise over the [2 x 81] matrices
  expect_true(all(totals[1L, ] == 10000L))
  expect_true(all(totals[2L, ] == 5000L))
})

test_that("run_sir flows equal the source-compartment deltas", {
  # Given the recorded trajectories and flows
  # When we reconstruct flows from the compartment differences
  # Then incidence[t] == S[t] - S[t+1] (the only outflow of S) and
  #      recovery[t]  == R[t+1] - R[t] (the only inflow to R), exactly.
  # Failure would mean the flow accounting is measuring the wrong step or window.
  nodes <- two_nodes()
  set.seed(4L)
  nticks <- 100L
  model <- run_sir(nodes, beta = 0.35, infectious_period = 8, nticks = nticks,
                   seed = c(15L, 3L))
  report <- attr(model, "report")

  s_drop <- report$S[, seq_len(nticks)] - report$S[, seq_len(nticks) + 1L]
  r_gain <- report$R[, seq_len(nticks) + 1L] - report$R[, seq_len(nticks)]
  expect_equal(report$incidence, s_drop)
  expect_equal(report$recovery, r_gain)
})

test_that("run_sir with zero seed produces no epidemic", {
  # Given a seed of 0 in every node
  # When the model runs
  # Then no agent is ever infected: incidence is zero throughout and S stays at
  #      the full population.
  # Failure would indicate spontaneous infection (force of infection computed
  #      from a non-zero phantom I, or mis-seeding).
  nodes <- two_nodes()
  set.seed(5L)
  model <- run_sir(nodes, beta = 0.5, infectious_period = 7, nticks = 30L,
                   seed = 0L)
  report <- attr(model, "report")

  expect_true(all(report$incidence == 0L))
  expect_equal(report$S[, ncol(report$S)], c(10000L, 5000L))
})

test_that("run_sir does not mix infection across unconnected nodes", {
  # Given infection seeded only in node 1 and no migration between nodes
  # When the model runs
  # Then node 2 (the unseeded node) never sees an infection.
  # Failure would mean the per-node force of infection is leaking across nodes.
  nodes <- two_nodes()
  set.seed(6L)
  model <- run_sir(nodes, beta = 0.6, infectious_period = 7, nticks = 60L,
                   seed = c(50L, 0L))
  report <- attr(model, "report")

  expect_true(all(report$I[2L, ] == 0L))
  expect_true(all(report$incidence[2L, ] == 0L))
})

test_that("run_sir attaches report and run metadata to the returned frame", {
  # Given a completed run
  # When we inspect the returned people frame's attributes
  # Then the report frame, model label, tick count, runtime, and parameters are
  #      all present and well-formed.
  # Failure would break the documented query surface of the return value.
  nodes <- two_nodes()
  set.seed(7L)
  model <- run_sir(nodes, beta = 0.3, infectious_period = 7, nticks = 25L,
                   seed = c(10L, 0L))

  expect_s3_class(attr(model, "report"), "LaserFrame")
  expect_equal(attr(model, "model"), "SIR")
  expect_equal(attr(model, "nticks"), 25L)
  expect_true(is.numeric(attr(model, "runtime")) && attr(model, "runtime") >= 0)
  expect_s3_class(attr(model, "parameters")$infectious_period, "Distribution")
  # Trajectory matrices are [n_nodes x (nticks + 1)]; flows are [n_nodes x nticks].
  expect_equal(dim(attr(model, "report")$I), c(2L, 26L))
  expect_equal(dim(attr(model, "report")$incidence), c(2L, 25L))
})

test_that("run_sir accepts a Distribution for the infectious period", {
  # Given a stochastic (gamma) infectious period rather than a fixed number
  # When the model runs
  # Then it completes and still conserves population — the Distribution path is
  #      wired through to the kernel and the seeded timers.
  # Failure would mean as_distribution() or the seed sampler rejected a real
  #      Distribution object.
  nodes <- two_nodes()
  set.seed(8L)
  model <- run_sir(nodes, beta = 0.4, infectious_period = dist_gamma(2, 3.5),
                   nticks = 50L, seed = c(20L, 0L))
  report <- attr(model, "report")

  expect_true(all(report$S + report$I + report$R == matrix(c(10000L, 5000L),
                                                            nrow = 2L,
                                                            ncol = 51L)))
})

test_that("run_sir validates its scenario table and arguments", {
  # Given malformed inputs
  # When run_sir is called
  # Then each contract violation raises an informative error rather than running
  #      a meaningless simulation.
  # Failure would let bad scenarios run silently and produce garbage.
  expect_error(run_sir(list(population = 100L), beta = 0.3,
                       infectious_period = 7, nticks = 10L),
               "data.frame")
  expect_error(run_sir(data.frame(pop = 100L), beta = 0.3,
                       infectious_period = 7, nticks = 10L),
               "population")
  expect_error(run_sir(two_nodes(), beta = 0.3, infectious_period = 7,
                       nticks = 0L),
               "nticks")
  expect_error(run_sir(two_nodes(), beta = 0.3, infectious_period = 7,
                       nticks = 10L, seed = c(1L, 2L, 3L)),
               "length 1 or")
  expect_error(run_sir(two_nodes(), beta = 0.3, infectious_period = 7,
                       nticks = 10L, seed = c(99999L, 0L)),
               "between 0 and")
})
