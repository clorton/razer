# Tests for the timed-transition step kernels. Each is a single pass that mutates the
# per-agent state/timer and RETURNS per-node transition counts (the caller applies them):
#   step_si   : M->S (waned), E->I (onset)
#   step_sir  : + I->absorbing (cleared)  [absorbing = S or R, parameterized]
#   step_sirs : + I->R (recovered, sets immunity timer), R->S (waned_r)
# Timers are u16. Written given-when-then.

states <- laser_states()                        # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]; R <- states[["R"]]; M <- states[["M"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }

test_that("step_si advances M->S and E->I on timer expiry", {
  # Given one node with M(t1), E(t1), and an untimed S and I
  # When step_si runs with a constant infectious period of 7
  # Then the M waner becomes S (timer 0), the E case becomes I with a fresh timer 7;
  #      S and I are untouched; the returned counts are waned=1, onset=1. Failure would
  #      mean a transition, the E->I timer draw, or the returned tally is wrong.
  state <- mk("u8", c(M, E, S, I)); timer <- mk("u16", c(1, 1, 0, 5)); nodeid <- mk("u16", c(0, 0, 0, 0))
  r <- step_si(state, timer, nodeid, 4L, 1L, dist_constant(7))

  expect_equal(state$values(), c(S, I, S, I))   # M->S, E->I; S,I untouched
  expect_equal(timer$values(), c(0, 7, 0, 5))   # E->I drew 7; I is terminal (timer untouched)
  expect_equal(r$waned, 1L); expect_equal(r$onset, 1L)
})

test_that("step_sir clears I to the parameterized absorbing state", {
  # Given an infectious agent expiring this tick
  # When step_sir runs with absorbing = S (SIS) vs absorbing = R (SIR)
  # Then the agent becomes S or R respectively, and `cleared` counts it. Failure would
  #      mean the absorbing-state parameter is ignored.
  to_s <- mk("u8", c(I)); step_sir(to_s, mk("u16", c(1)), mk("u16", c(0)), 1L, 1L, dist_constant(5), S)
  expect_equal(to_s$values(), S)

  to_r <- mk("u8", c(I)); r <- step_sir(to_r, mk("u16", c(1)), mk("u16", c(0)), 1L, 1L, dist_constant(5), R)
  expect_equal(to_r$values(), R)
  expect_equal(r$cleared, 1L); expect_equal(r$onset, 0L); expect_equal(r$waned, 0L)
})

test_that("step_sirs does M->S, E->I, I->R (with waning), and R->S", {
  # Given one of each timed state, all expiring this tick (M, E, I, R), inf period 6,
  #       immunity period 30
  # When step_sirs runs
  # Then M->S, E->I (timer 6), I->R (timer 30, the immunity clock), R->S (timer 0); the
  #      returned counts are waned=onset=recovered=waned_r=1. Failure would mean a
  #      transition, the immunity-timer set on I->R, or a tally is wrong.
  state <- mk("u8", c(M, E, I, R)); timer <- mk("u16", c(1, 1, 1, 1)); nodeid <- mk("u16", c(0, 0, 0, 0))
  r <- step_sirs(state, timer, nodeid, 4L, 1L, dist_constant(6), dist_constant(30))

  expect_equal(state$values(), c(S, I, R, S))
  expect_equal(timer$values(), c(0, 6, 30, 0))  # E->I draws 6; I->R sets the immunity 30
  expect_equal(r$waned, 1L); expect_equal(r$onset, 1L)
  expect_equal(r$recovered, 1L); expect_equal(r$waned_r, 1L)
})

test_that("step kernels only decrement timers that have not yet expired", {
  # Given timed agents whose timers do not expire this tick
  # When step_sir runs
  # Then no state changes, each timer is decremented by 1, and all counts are 0.
  state <- mk("u8", c(M, E, I)); timer <- mk("u16", c(3, 2, 5)); nodeid <- mk("u16", c(0, 0, 0))
  r <- step_sir(state, timer, nodeid, 3L, 1L, dist_constant(5), R)

  expect_equal(state$values(), c(M, E, I))
  expect_equal(timer$values(), c(2, 1, 4))
  expect_equal(r$waned, 0L); expect_equal(r$onset, 0L); expect_equal(r$cleared, 0L)
})

test_that("step_sir's per-node counts match a census at scale", {
  # Given 900,000 agents (300k each M/E/I, all timer 1) over 3 nodes
  # When step_sir(absorbing = R) runs (work split across cores)
  # Then every M->S, E->I, I->R fires; the returned per-node counts equal serial
  #      tabulates of the resulting agent states. Failure would expose a reduction race.
  set.seed(9L)
  per <- 300000L; nn <- 3L
  nid <- sample.int(nn, 3L * per, replace = TRUE) - 1L
  st  <- c(rep(M, per), rep(E, per), rep(I, per))
  state <- mk("u8", st); timer <- mk("u16", rep(1L, 3L * per)); nodeid <- mk("u16", nid)
  r <- step_sir(state, timer, nodeid, 3L * per, nn, dist_constant(6), R)

  expect_equal(r$waned,   as.integer(tabulate(nid[st == M] + 1L, nn)))   # all M waned
  expect_equal(r$onset,   as.integer(tabulate(nid[st == E] + 1L, nn)))   # all E onset
  expect_equal(r$cleared, as.integer(tabulate(nid[st == I] + 1L, nn)))   # all I cleared
  expect_true(all(state$values()[st == M] == S))
  expect_true(all(state$values()[st == E] == I))
  expect_true(all(state$values()[st == I] == R))
})
