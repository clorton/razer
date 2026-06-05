# Tests for transmission(): the per-tick infection kernel. It converts the per-node
# FOI rate to a per-node probability (1 - exp(-foi), once per node), moves
# susceptibles into a receiving state (`to_state` = I for SIR, E for SEIR), sets
# their timer from a Distribution, and applies the S->to_state delta to census
# column tick+1 (plus the per-node incidence flow). Assumes the census was already
# carried forward this tick. Written given-when-then. The parallel-accumulation test
# checks the per-core tally against a serial census of the resulting agent states.

states <- laser_states()
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }

test_that("transmission into I (SIR) applies the S->I delta and flow", {
  # Given 5 susceptibles in 2 nodes, a census carried to tick 1, and a FOI that is
  #       huge in node 0 (p ~ 1) and zero in node 1, with to_state = I
  # When transmission runs for tick 0
  # Then node-0 susceptibles become I (with a drawn infectious timer), node 1 stays
  #       S; S/I census at tick 1 and incidence reflect the new infections; and the
  #       census equals a direct census of the agents.
  # Failure would mean the probability, the receiving state, the delta, or the flow
  #       is wrong.
  state  <- mk("u8",  rep(S, 5))
  timer  <- mk("u8",  rep(0, 5))
  nodeid <- mk("u16", c(0, 0, 1, 1, 1))
  Sc <- allocate_vector("i32", 3L, 2L); Sc$set(c(2, 3, 2, 3, 0, 0))   # col0 & col1 (carried)
  Ic <- allocate_vector("i32", 3L, 2L); Ic$set(rep(0L, 6L))
  foi <- allocate_vector("f64", 2L, 2L); foi$set(c(100, 0, 0, 0))
  incidence <- allocate_vector("i32", 2L, 2L)

  transmission(state, timer, nodeid, 5L, foi, Sc, Ic, incidence, 0L, I, dist_constant(7))

  sv <- state$values()
  expect_equal(Sc$values()[2L, ], c(0, 3))
  expect_equal(Ic$values()[2L, ], c(2, 0))
  expect_equal(incidence$values()[1L, ], c(2, 0))
  expect_true(all(timer$values()[sv == I] == 7))
  expect_equal(Ic$values()[2L, ], as.numeric(tabulate(nodeid$values()[sv == I] + 1L, 2L)))
})

test_that("transmission into E (SEIR) moves susceptibles to Exposed with an incubation timer", {
  # Given the same setup but to_state = E and an E census instead of I
  # When transmission runs
  # Then node-0 susceptibles become E (not I), the E census and incidence rise, the
  #       infectious census is untouched, and the timer holds the incubation draw.
  # Failure would mean the receiving state or its census is hard-wired to I.
  state  <- mk("u8",  rep(S, 5))
  timer  <- mk("u8",  rep(0, 5))
  nodeid <- mk("u16", c(0, 0, 1, 1, 1))
  Sc <- allocate_vector("i32", 3L, 2L); Sc$set(c(2, 3, 2, 3, 0, 0))
  Ec <- allocate_vector("i32", 3L, 2L); Ec$set(rep(0L, 6L))
  Ic <- allocate_vector("i32", 3L, 2L); Ic$set(rep(0L, 6L))   # should stay zero
  foi <- allocate_vector("f64", 2L, 2L); foi$set(c(100, 0, 0, 0))
  incidence <- allocate_vector("i32", 2L, 2L)

  transmission(state, timer, nodeid, 5L, foi, Sc, Ec, incidence, 0L, E, dist_constant(4))

  sv <- state$values()
  expect_equal(sum(sv == E), 2L)                     # two moved to Exposed
  expect_equal(sum(sv == I), 0L)                     # none went straight to I
  expect_equal(Sc$values()[2L, ], c(0, 3))
  expect_equal(Ec$values()[2L, ], c(2, 0))           # E census rose
  expect_equal(Ic$values()[2L, ], c(0, 0))           # I census untouched
  expect_true(all(timer$values()[sv == E] == 4))     # incubation timer set
})

test_that("transmission's parallel tally matches a census at scale", {
  # Given 1,000,000 susceptibles over 200 nodes, census carried to tick 1, and a
  #       uniform FOI giving p ~ 0.5 per node, to_state = I
  # When transmission runs (work split across cores with per-core node buffers)
  # Then — regardless of the stochastic outcome — the incidence total equals the
  #       number now infectious, the S/I census equals a serial tabulate of the
  #       agents by node, and S + I is conserved per node.
  # Failure would expose a race or reduction bug in the parallel accumulation.
  set.seed(2L)
  n_agents <- 1000000L; n_nodes <- 200L
  nid <- sample.int(n_nodes, n_agents, replace = TRUE) - 1L
  state  <- mk("u8",  rep(S, n_agents))
  timer  <- mk("u8",  rep(0, n_agents))
  nodeid <- mk("u16", nid)
  per_node <- as.numeric(tabulate(nid + 1L, n_nodes))
  Sc <- allocate_vector("i32", 2L, n_nodes); Sc$set(c(per_node, per_node))
  Ic <- allocate_vector("i32", 2L, n_nodes); Ic$set(rep(0L, 2L * n_nodes))
  foi <- allocate_vector("f64", 1L, n_nodes); foi$set(rep(0.7, n_nodes))
  incidence <- allocate_vector("i32", 1L, n_nodes)

  transmission(state, timer, nodeid, n_agents, foi, Sc, Ic, incidence, 0L, I, dist_constant(5))

  sv <- state$values(); n_inf <- sum(sv == I)
  expect_equal(sum(incidence$values()), n_inf)
  expect_equal(Ic$values()[2L, ], as.numeric(tabulate(nodeid$values()[sv == I] + 1L, n_nodes)))
  expect_equal(Sc$values()[2L, ], as.numeric(tabulate(nodeid$values()[sv == S] + 1L, n_nodes)))
  expect_equal(Sc$values()[2L, ] + Ic$values()[2L, ], per_node)
  expect_gt(n_inf, 0); expect_lt(n_inf, n_agents)
})

test_that("transmission validates the tick range and to_state", {
  # Given a tick beyond the FOI buffer and an out-of-range state code
  # When transmission is called
  # Then it errors rather than reading out of bounds or casting a bad state code.
  # Failure would risk out-of-bounds access.
  state <- mk("u8", c(S, S)); timer <- mk("u8", c(0, 0)); nodeid <- mk("u16", c(0, 1))
  Sc <- allocate_vector("i32", 3L, 2L); Ic <- allocate_vector("i32", 3L, 2L)
  foi <- allocate_vector("f64", 2L, 2L); incidence <- allocate_vector("i32", 2L, 2L)
  expect_error(transmission(state, timer, nodeid, 2L, foi, Sc, Ic, incidence, 2L, I, dist_constant(7)),
               "out of range")
  expect_error(transmission(state, timer, nodeid, 2L, foi, Sc, Ic, incidence, 0L, 300L, dist_constant(7)),
               "state code")
})
