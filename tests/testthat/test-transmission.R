# Tests for the transmission kernels. Both convert column `tick` of `foi` to a per-node
# probability 1-exp(-foi), move susceptibles into the disease chain, and RETURN the
# per-node count of new infections (the caller applies the census deltas):
#   transmission     — S->to_state (E or I), sets a u16 timer from a duration.
#   transmission_si  — S->I, I absorbing (no timer) — the SI model.
# The RNG is thread-local (not R-seedable), so checks are deterministic-at-the-extremes
# (p~1 / p=0) or statistical at scale. Written given-when-then.

states <- laser_states()                        # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }
foi1 <- function(v) { g <- allocate_vector("f64", 1L, length(v)); g$set(v); g }  # one-tick foi

test_that("transmission moves S->to_state, sets a u16 timer, and returns counts", {
  # Given 4 susceptibles in one node and a huge FOI (p ~ 1)
  # When transmission runs with to_state = E and a constant incubation of 9
  # Then all 4 become E with timer 9, and the returned per-node count is 4. Failure would
  #      mean the move, the timer, or the returned tally is wrong.
  state <- mk("u8", rep(S, 4L)); timer <- mk("u16", rep(0L, 4L)); nodeid <- mk("u16", rep(0L, 4L))
  inf <- transmission(state, timer, nodeid, 4L, foi1(100), 0L, E, dist_constant(9))

  expect_equal(inf, 4L)                          # per-node infections (1 node)
  expect_true(all(state$values() == E))
  expect_true(all(timer$values() == 9L))
})

test_that("transmission with zero FOI infects nobody and returns zero", {
  # Given susceptibles and a zero FOI
  # When transmission runs
  # Then no agent changes and the returned count is 0.
  state <- mk("u8", rep(S, 3L)); timer <- mk("u16", rep(0L, 3L)); nodeid <- mk("u16", rep(0L, 3L))
  inf <- transmission(state, timer, nodeid, 3L, foi1(0), 0L, I, dist_constant(9))

  expect_equal(inf, 0L)
  expect_true(all(state$values() == S))
})

test_that("transmission_si moves S->I absorbing with no timer, and returns counts", {
  # Given susceptibles and p ~ 1
  # When transmission_si runs
  # Then all become I and the count is returned (no timer is set — I is terminal in SI).
  state <- mk("u8", rep(S, 5L)); nodeid <- mk("u16", rep(0L, 5L))
  inf <- transmission_si(state, nodeid, 5L, foi1(100), 0L)

  expect_equal(inf, 5L)
  expect_true(all(state$values() == I))
})

test_that("transmission's per-node counts match the agent state changes at scale", {
  # Given 1,000,000 susceptibles over 100 nodes and a moderate per-node FOI (p ~ 0.5)
  # When transmission runs (work split across cores)
  # Then the returned per-node counts equal a serial tabulate of the now-E agents, and
  #      the total is strictly between 0 and n. Failure would expose a race in the tally.
  set.seed(5L)
  n_agents <- 1000000L; n_nodes <- 100L
  nid <- sample.int(n_nodes, n_agents, replace = TRUE) - 1L
  state <- mk("u8", rep(S, n_agents)); timer <- mk("u16", rep(0L, n_agents)); nodeid <- mk("u16", nid)
  inf <- transmission(state, timer, nodeid, n_agents, foi1(rep(log(2), n_nodes)), 0L, E, dist_constant(9))

  exposed <- as.integer(tabulate(nodeid$values()[state$values() == E] + 1L, n_nodes))
  expect_equal(inf, exposed)
  expect_gt(sum(inf), 0); expect_lt(sum(inf), n_agents)
})
