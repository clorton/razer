# Generalized timer-expiry kernels: step_timer_expire (absorbing destination)
# and step_timer_expire_set (timed destination).
#
# These mirror laser-generic's `nb_timer_update` and `nb_timer_update_timer_set`
# and are the engines behind the four named kernels:
#   step_infectious_is(people)            == step_timer_expire(people, I, S)
#   step_recovered_rs(people)             == step_timer_expire(people, R, S)
#   step_exposed_ei(people, d)            == step_timer_expire_set(people, E, I, d)
#   step_infectious_ir(people, d)         == step_timer_expire_set(people, I, R, d)
#
# State codes come from laser_states(): S=0, E=1, I=2, R=3, D=-1.
#
# testthat idioms: `test_that("desc", { ... })` blocks; `expect_equal` compares
# values; `L` suffixes are integer literals.

# `[["name"]]` extracts one element of the named vector by name (no partial match),
# binding module-level constants S/E/I/R used throughout the assertions.
S <- laser_states()[["S"]]
E <- laser_states()[["E"]]
I <- laser_states()[["I"]]
R <- laser_states()[["R"]]

# Build a small people frame with a given vector of states and timers.
make_people <- function(states, timers) {
  n <- length(states)
  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)
  ppl$state <- as.integer(states)   # `$prop <- v` writes the whole column; as.integer coerces type
  ppl$timer <- as.integer(timers)
  ppl   # last expression is the return value
}

# ── Absorbing-state transition: step_timer_expire ───────────────────────────────

test_that("step_timer_expire decrements and transitions only the from_state on expiry", {
  # Given: three I agents with timers 1, 2, 3 and one S bystander with timer 5
  # When:  one application of step_timer_expire(I -> S)
  # Then:  only the timer-1 I agent transitions to S (timer reset to 0); the
  #        others decrement but stay I; the S bystander is untouched.
  # Failure would mean the from_state guard or the expiry test is wrong.
  ppl <- make_people(c(I, I, I, S), c(1L, 2L, 3L, 5L))
  step_timer_expire(ppl, from_state = I, to_state = S)
  expect_equal(ppl$state, c(S, I, I, S))
  expect_equal(ppl$timer, c(0L, 1L, 2L, 5L))  # S bystander's timer never decremented
})

test_that("step_timer_expire leaves an absorbing destination with timer 0", {
  # Given: an I agent that expires into S
  # When:  the transition fires
  # Then:  the destination carries no duration of its own (timer == 0), so a
  #        subsequent decrement does not immediately re-fire.
  ppl <- make_people(I, 1L)
  step_timer_expire(ppl, from_state = I, to_state = S)
  expect_equal(ppl$state, S)
  expect_equal(ppl$timer, 0L)
})

# ── Timed-destination transition: step_timer_expire_set ─────────────────────────

test_that("step_timer_expire_set assigns a fresh destination timer on expiry", {
  # Given: an E agent with timer 1 transitioning to I with a constant 7-tick
  #        infectious period
  # When:  step_timer_expire_set(E -> I, dist_constant(7)) fires
  # Then:  the agent is now I with a fresh timer of 7 (the destination's duration)
  # Failure would mean the new-timer draw is not applied at transition.
  ppl <- make_people(E, 1L)
  step_timer_expire_set(ppl, from_state = E, to_state = I,
                        duration_dist = dist_constant(7))
  expect_equal(ppl$state, I)
  expect_equal(ppl$timer, 7L)
})

test_that("step_timer_expire_set clamps a sub-tick draw to a minimum of one tick", {
  # Given: an I agent expiring into R with a degenerate 0-tick duration
  # When:  the transition fires
  # Then:  the destination timer is clamped to 1 (a timed state lasts >= 1 tick)
  ppl <- make_people(I, 1L)
  step_timer_expire_set(ppl, from_state = I, to_state = R,
                        duration_dist = dist_constant(0))
  expect_equal(ppl$state, R)
  expect_equal(ppl$timer, 1L)   # max(round(0), 1)
})

# ── Equivalence of the named wrappers to the generalized kernels ────────────────

test_that("step_infectious_is equals step_timer_expire(I, S)", {
  # Given: identical I-with-timer frames driven by the named and generalized kernels
  # When:  each frame is advanced one tick
  # Then:  state and timer vectors are identical, proving the wrapper delegates
  #        to the same core loop.
  states <- c(I, I, I, S)
  timers <- c(1L, 2L, 1L, 4L)
  a <- make_people(states, timers); step_infectious_is(a)
  b <- make_people(states, timers); step_timer_expire(b, from_state = I, to_state = S)
  expect_equal(a$state, b$state)
  expect_equal(a$timer, b$timer)
})

test_that("step_recovered_rs equals step_timer_expire(R, S)", {
  states <- c(R, R, S, I)
  timers <- c(1L, 3L, 0L, 2L)
  a <- make_people(states, timers); step_recovered_rs(a)
  b <- make_people(states, timers); step_timer_expire(b, from_state = R, to_state = S)
  expect_equal(a$state, b$state)
  expect_equal(a$timer, b$timer)
})

test_that("step_exposed_ei equals step_timer_expire_set(E, I, dist) under a constant duration", {
  # A constant distribution makes the draw deterministic, so the two kernels must
  # agree exactly (no RNG divergence to worry about).
  states <- c(E, E, I, S)
  timers <- c(1L, 2L, 5L, 0L)
  a <- make_people(states, timers); step_exposed_ei(a, inf_dist = dist_constant(9))
  b <- make_people(states, timers)
  step_timer_expire_set(b, from_state = E, to_state = I, duration_dist = dist_constant(9))
  expect_equal(a$state, b$state)
  expect_equal(a$timer, b$timer)
})

test_that("step_infectious_ir equals step_timer_expire_set(I, R, dist) under a constant duration", {
  states <- c(I, I, R, S)
  timers <- c(1L, 4L, 2L, 0L)
  a <- make_people(states, timers); step_infectious_ir(a, imm_dist = dist_constant(12))
  b <- make_people(states, timers)
  step_timer_expire_set(b, from_state = I, to_state = R, duration_dist = dist_constant(12))
  expect_equal(a$state, b$state)
  expect_equal(a$timer, b$timer)
})

# ── Composability: a full SEIRS chain built only from generalized kernels ───────

test_that("a single E agent walks E -> I -> R -> S using only generalized kernels", {
  # Given: one E agent, exp = 2, inf = 3, imm = 2 ticks
  # When:  applied downstream-first each tick (R->S, I->R, E->I)
  # Then:  it occupies E for 2 ticks, I for 3, R for 2, then returns to S,
  #        demonstrating the two generalized kernels compose into a full model.
  ppl <- make_people(E, 2L)
  observed <- integer(0)   # an empty integer vector, used as a growable accumulator
  for (tick in seq_len(8L)) {
    step_timer_expire(ppl, from_state = R, to_state = S)               # R -> S (absorbing)
    step_timer_expire_set(ppl, from_state = I, to_state = R,
                          duration_dist = dist_constant(2))            # I -> R (timed)
    step_timer_expire_set(ppl, from_state = E, to_state = I,
                          duration_dist = dist_constant(3))            # E -> I (timed)
    observed <- c(observed, ppl$state[1L])   # c() appends by building a new vector each tick
  }
  # Ticks:        1  2  3  4  5  6  7  8
  # Expected:     E  I  I  I  R  R  S  S
  expect_equal(observed, c(E, I, I, I, R, R, S, S))
})
