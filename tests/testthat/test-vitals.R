# Tests for constant_pop_vitals_sir(): constant-population SIR vital dynamics. Each
# agent dies with probability (1 - exp(-daily_rate)) and is immediately reborn
# susceptible (state -> S, timer -> 0); every event is recorded as both a birth and a
# death (equal under constant population), and the S/I/R census is kept in sync
# (deaths out of I/R move those counts to S; the census delta lands at column tick+1).
# Written given-when-then. The parallel test checks the per-core tally against a
# serial census of the resulting agent states.

states <- laser_states()
S <- states[["S"]]; I <- states[["I"]]; R <- states[["R"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }
# A 2-column i32 census with column 0 AND column 1 = `col0` (column 1 simulates the
# caller's carry_forward for tick 0).
carried <- function(col0) {
  buf <- allocate_vector("i32", 2L, length(col0)); buf$set(c(col0, col0)); buf
}

test_that("constant_pop_vitals_sir reverts the dead to S and keeps the census in sync", {
  # Given node 0 holding I, I, R (S=0, I=2, R=1) with a huge death rate (p ~ 1), and
  #       node 1 holding S, R (S=1, R=1) with a zero death rate, census carried to tick 1
  # When constant_pop_vitals_sir runs for tick 0
  # Then every node-0 agent dies and is reborn susceptible; node-0 census becomes
  #       S=3 (the two ex-I and one ex-R), I=0, R=0; node 1 is untouched; births ==
  #       deaths == per-node events; S+I+R is conserved and matches the agents.
  # Failure would mean the reversion, the per-compartment census delta, or the tally
  #       is wrong.
  state  <- mk("u8",  c(I, I, R, S, R))
  timer  <- mk("u16",  c(7, 3, 0, 0, 5))
  nodeid <- mk("u16", c(0, 0, 0, 1, 1))
  rate   <- allocate_vector("f64", 2L, 2L); rate$set(c(100, 0, 0, 0))
  Sc <- carried(c(0, 1)); Ic <- carried(c(2, 0)); Rc <- carried(c(1, 1))
  births <- allocate_vector("i32", 1L, 2L); deaths <- allocate_vector("i32", 1L, 2L)

  constant_pop_vitals_sir(state, timer, nodeid, 5L, rate, Sc, Ic, Rc, births, deaths, 0L)

  sv <- state$values()
  expect_equal(state$values(), c(S, S, S, S, R))   # node0 all S; node1 unchanged
  expect_equal(timer$values(), c(0, 0, 0, 0, 5))
  expect_equal(Sc$values()[2L, ], c(3, 1))         # node0: 0 + (2 ex-I + 1 ex-R)
  expect_equal(Ic$values()[2L, ], c(0, 0))
  expect_equal(Rc$values()[2L, ], c(0, 1))
  expect_equal(births$values(), c(3, 0))
  expect_equal(deaths$values(), c(3, 0))
  # S+I+R conserved per node, and equal to a direct agent census
  expect_equal(Sc$values()[2L, ] + Ic$values()[2L, ] + Rc$values()[2L, ], c(3, 2))
  expect_equal(Sc$values()[2L, 1], sum(sv[nodeid$values() == 0] == S))
  expect_equal(Rc$values()[2L, 1], sum(sv[nodeid$values() == 0] == R))
})

test_that("constant_pop_vitals_sir does nothing at a zero death rate", {
  # Given a zero death rate everywhere, census carried to tick 1
  # When constant_pop_vitals_sir runs
  # Then no agent changes, the census is unchanged, and births/deaths are zero.
  # Failure would mean spurious deaths or census drift at zero hazard.
  state  <- mk("u8",  c(I, R, S)); timer <- mk("u16", c(4, 0, 0)); nodeid <- mk("u16", c(0, 0, 1))
  rate   <- allocate_vector("f64", 2L, 2L)
  Sc <- carried(c(1, 1)); Ic <- carried(c(1, 0)); Rc <- carried(c(1, 0))
  births <- allocate_vector("i32", 1L, 2L); deaths <- allocate_vector("i32", 1L, 2L)

  constant_pop_vitals_sir(state, timer, nodeid, 3L, rate, Sc, Ic, Rc, births, deaths, 0L)

  expect_equal(state$values(), c(I, R, S))
  expect_equal(Sc$values()[2L, ], c(1, 1)); expect_equal(Ic$values()[2L, ], c(1, 0))
  expect_equal(Rc$values()[2L, ], c(1, 0))
  expect_equal(births$values(), c(0, 0)); expect_equal(deaths$values(), c(0, 0))
})

test_that("constant_pop_vitals_sir's parallel tally and census match a census at scale", {
  # Given 1,000,000 infectious agents over 200 nodes (census I = per-node, S = R = 0,
  #       carried to tick 1) and a moderate death rate (p ~ 0.5)
  # When constant_pop_vitals_sir runs (work split across cores with per-core buffers)
  # Then every death is an I->S move, so I drops and S rises by the deaths; births ==
  #       deaths == the per-node reborn count == a serial tabulate of the now-S agents;
  #       the census equals a direct agent census and S+I is conserved.
  # Failure would expose a race or reduction bug in the parallel accumulation.
  set.seed(3L)
  n_agents <- 1000000L; n_nodes <- 200L
  nid <- sample.int(n_nodes, n_agents, replace = TRUE) - 1L
  state  <- mk("u8",  rep(I, n_agents)); timer <- mk("u16", rep(7L, n_agents)); nodeid <- mk("u16", nid)
  per_node <- as.numeric(tabulate(nid + 1L, n_nodes))
  rate   <- allocate_vector("f64", 1L, n_nodes); rate$set(rep(0.7, n_nodes))
  Sc <- carried(rep(0, n_nodes)); Ic <- carried(per_node); Rc <- carried(rep(0, n_nodes))
  births <- allocate_vector("i32", 1L, n_nodes); deaths <- allocate_vector("i32", 1L, n_nodes)

  constant_pop_vitals_sir(state, timer, nodeid, n_agents, rate, Sc, Ic, Rc, births, deaths, 0L)

  sv <- state$values()
  reborn <- as.numeric(tabulate(nodeid$values()[sv == S] + 1L, n_nodes))
  expect_equal(births$values(), reborn)                          # events == now-S agents
  expect_equal(deaths$values(), births$values())                 # births == deaths
  expect_equal(Sc$values()[2L, ], reborn)                        # S census == now-S agents
  expect_equal(Ic$values()[2L, ], per_node - reborn)             # I census == survivors
  expect_equal(Sc$values()[2L, ] + Ic$values()[2L, ], per_node)  # S + I conserved
  n_reborn <- sum(sv == S)
  expect_gt(n_reborn, 0); expect_lt(n_reborn, n_agents)
})

test_that("constant_pop_vitals_sir validates the tick range and shapes", {
  # Given a tick whose census column tick+1 does not exist
  # When constant_pop_vitals_sir is called
  # Then it errors rather than writing out of bounds.
  # Failure would risk out-of-bounds access.
  state <- mk("u8", c(I, I)); timer <- mk("u16", c(1, 1)); nodeid <- mk("u16", c(0, 1))
  rate <- allocate_vector("f64", 2L, 2L)
  Sc <- carried(c(1, 1)); Ic <- carried(c(1, 1)); Rc <- carried(c(0, 0))
  births <- allocate_vector("i32", 2L, 2L); deaths <- allocate_vector("i32", 2L, 2L)
  expect_error(
    constant_pop_vitals_sir(state, timer, nodeid, 2L, rate, Sc, Ic, Rc, births, deaths, 1L),
    "out of range")
})
