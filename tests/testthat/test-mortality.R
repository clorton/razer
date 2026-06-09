# Tests for mortality(): natural-mortality step that retires agents whose date of death
# `dod` (an absolute tick) has arrived. It mutates the per-agent state (D = 255) and
# RETURNS per-node death counts broken down by source state, `list(m, s, e, i, r)`;
# the caller decrements those census states. Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]; R <- states[["R"]]; M <- states[["M"]]
D_U8 <- 255L                                   # how D (-1) reads back from a u8 column

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }

test_that("mortality retires due agents and returns deaths by state", {
  # Given (tick 0, two nodes) node 0 = {I(dod0), S(dod0), R(dod5)}, node 1 = {M(dod0),
  #       S(dod3)}
  # When mortality runs for tick 0
  # Then the agents with dod <= 0 (the I and S in node 0, the M in node 1) become D; the
  #       R (dod 5) and node-1 S (dod 3) survive; the returned counts are i=(1,0),
  #       s=(1,0), m=(0,1), with e and r zero. Failure would mean the due-date test or the
  #       per-state tally is wrong.
  state  <- mk("u8",  c(I, S, R, M, S))
  dod    <- mk("u32", c(0, 0, 5, 0, 3))
  nodeid <- mk("u16", c(0, 0, 0, 1, 1))

  d <- mortality(state, dod, nodeid, 5L, 2L, 0L)

  expect_equal(state$values(), c(D_U8, D_U8, R, D_U8, S))
  expect_equal(d$i, c(1L, 0L)); expect_equal(d$s, c(1L, 0L)); expect_equal(d$m, c(0L, 1L))
  expect_equal(d$e, c(0L, 0L)); expect_equal(d$r, c(0L, 0L))
})

test_that("mortality leaves agents whose dod is in the future untouched", {
  # Given agents all with dod > tick
  # When mortality runs
  # Then nobody dies and every returned count is zero.
  state  <- mk("u8", c(S, I, R)); dod <- mk("u32", c(5, 9, 7)); nodeid <- mk("u16", c(0, 0, 0))
  d <- mortality(state, dod, nodeid, 3L, 1L, 0L)

  expect_equal(state$values(), c(S, I, R))
  expect_equal(d$s, 0L); expect_equal(d$i, 0L); expect_equal(d$r, 0L)
})

test_that("mortality does not re-count already-deceased agents", {
  # Given an agent already D (255) alongside a due living agent
  # When mortality runs
  # Then only the living agent is retired and counted (no double-count).
  state  <- mk("u8", c(D_U8, I)); dod <- mk("u32", c(0, 0)); nodeid <- mk("u16", c(0, 0))
  d <- mortality(state, dod, nodeid, 2L, 1L, 0L)

  expect_equal(state$values(), c(D_U8, D_U8))
  expect_equal(d$i, 1L)                          # the one live I retired
})

test_that("mortality's per-node counts match a census at scale", {
  # Given 1,000,000 infectious agents over 50 nodes, half with dod = 0 (due) and half in
  #       the future
  # When mortality runs for tick 0 (work split across cores)
  # Then the returned per-node I-death counts equal a serial tabulate of the now-D agents.
  set.seed(13L)
  n_agents <- 1000000L; n_nodes <- 50L
  nid <- sample.int(n_nodes, n_agents, replace = TRUE) - 1L
  due <- sample(c(0L, 10L), n_agents, replace = TRUE)   # dod 0 (die) or 10 (survive tick 0)
  state <- mk("u8", rep(I, n_agents)); dod <- mk("u32", due); nodeid <- mk("u16", nid)

  d <- mortality(state, dod, nodeid, n_agents, n_nodes, 0L)

  dead <- as.integer(tabulate(nid[state$values() == D_U8] + 1L, n_nodes))
  expect_equal(d$i, dead)
  expect_gt(sum(d$i), 0); expect_lt(sum(d$i), n_agents)
})
