# Tests for transmission_u16(): the measles S->E transmission step. Identical to
# transmission() but writes a uint16 timer (the incubation clock). It converts foi[tick]
# to a per-node probability 1-exp(-foi), moves susceptibles into the receiving state
# (E), draws their u16 timer from `duration`, and applies the S-down / E-up census delta
# at column tick+1 (recording incidence at column tick). Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }
carried <- function(col0) {
  buf <- allocate_vector("i32", 2L, length(col0)); buf$set(c(col0, col0)); buf
}

test_that("transmission_u16 with a huge FOI exposes every susceptible (S->E)", {
  # Given 4 susceptibles in node 0 and a very large FOI (p ~ 1), census carried to tick 1
  # When transmission_u16 runs for tick 0 with a constant incubation of 9 and to_state = E
  # Then every S becomes E with a uint16 timer of 9; the census moves S 4->0, E 0->4; and
  #      incidence records 4. Failure would mean S->E, the u16 timer, or the delta is wrong.
  state  <- mk("u8",  rep(S, 4L)); timer <- mk("u16", rep(0L, 4L)); nodeid <- mk("u16", rep(0L, 4L))
  foi    <- allocate_vector("f64", 1L, 1L); foi$set(100)        # p = 1 - exp(-100) ~ 1
  Sc <- carried(4); Ec <- carried(0); incidence <- allocate_vector("i32", 1L, 1L)

  transmission_u16(state, timer, nodeid, 4L, foi, Sc, Ec, incidence, 0L, E, dist_constant(9))

  expect_true(all(state$values() == E))
  expect_true(all(timer$values() == 9L))                       # u16 incubation timer
  expect_equal(Sc$values()[2L, ], 0); expect_equal(Ec$values()[2L, ], 4)
  expect_equal(incidence$values(), 4)
})

test_that("transmission_u16 with zero FOI exposes nobody", {
  # Given susceptibles and a zero FOI
  # When transmission_u16 runs
  # Then no agent changes and the census/incidence stay put. Failure would mean spurious
  #      exposures at zero force of infection.
  state  <- mk("u8",  rep(S, 3L)); timer <- mk("u16", rep(0L, 3L)); nodeid <- mk("u16", rep(0L, 3L))
  foi    <- allocate_vector("f64", 1L, 1L)                     # 0
  Sc <- carried(3); Ec <- carried(0); incidence <- allocate_vector("i32", 1L, 1L)

  transmission_u16(state, timer, nodeid, 3L, foi, Sc, Ec, incidence, 0L, E, dist_constant(9))

  expect_true(all(state$values() == S))
  expect_equal(Sc$values()[2L, ], 3); expect_equal(Ec$values()[2L, ], 0)
  expect_equal(incidence$values(), 0)
})

test_that("transmission_u16's per-node exposures match a census at scale", {
  # Given 1,000,000 susceptibles over 100 nodes and a moderate per-node FOI (p ~ 0.5),
  #       census carried to tick 1
  # When transmission_u16 runs (work split across cores)
  # Then the per-node E census equals a serial tabulate of the now-E agents, S falls by
  #      the same, and incidence matches; total S+E is conserved. Failure would expose a
  #      race in the parallel accumulation.
  set.seed(11L)
  n_agents <- 1000000L; n_nodes <- 100L
  nid <- sample.int(n_nodes, n_agents, replace = TRUE) - 1L
  state  <- mk("u8",  rep(S, n_agents)); timer <- mk("u16", rep(0L, n_agents)); nodeid <- mk("u16", nid)
  foi    <- allocate_vector("f64", 1L, n_nodes); foi$set(rep(log(2), n_nodes))   # p = 0.5
  Sc <- carried(as.numeric(tabulate(nid + 1L, n_nodes))); Ec <- carried(rep(0, n_nodes))
  incidence <- allocate_vector("i32", 1L, n_nodes)

  transmission_u16(state, timer, nodeid, n_agents, foi, Sc, Ec, incidence, 0L, E, dist_constant(9))

  sv <- state$values()
  exposed <- as.numeric(tabulate(nodeid$values()[sv == E] + 1L, n_nodes))
  expect_equal(Ec$values()[2L, ], exposed)
  expect_equal(incidence$values(), exposed)
  expect_equal(Sc$values()[2L, ] + Ec$values()[2L, ], as.numeric(tabulate(nid + 1L, n_nodes)))
  n_exp <- sum(sv == E); expect_gt(n_exp, 0); expect_lt(n_exp, n_agents)
})
