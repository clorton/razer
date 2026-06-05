# Tests for sir_step(): the per-tick SIR recovery kernel over Column buffers. It
# decrements the timer of each infectious agent and, when the timer hits zero,
# moves the agent to Recovered and tallies a recovery for its node in the current
# tick's column of a 2-D report. Written given-when-then; `L` marks an integer
# literal. Deterministic, so exact assertions.

states <- laser_states()                 # c(S=0, E=1, I=2, R=3, D=-1)
S <- states[["S"]]; I <- states[["I"]]; R <- states[["R"]]

make <- function(dtype, x) {
  col <- allocate_scalar(dtype, length(x))
  col$set(x)
  col
}

test_that("sir_step recovers expired infectious agents and tallies them by node", {
  # Given 4 agents in 2 nodes — three infectious (timers 1, 2, 1) and one
  #       susceptible — and an empty 2-node x 3-tick recoveries report
  # When one sir_step runs for tick 0
  # Then the two agents whose timer reaches 0 become R and are counted in their
  #      nodes' tick-0 column; the agent with timer 2 stays I (timer now 1); the
  #      susceptible is untouched; other tick columns stay zero.
  # Failure would mean the recovery transition, the timer countdown, or the
  #      per-node/per-tick tally is wrong.
  state  <- make("u8",  c(I, I, S, I))
  timer  <- make("u8",  c(1, 2, 0, 1))
  nodeid <- make("u16", c(0, 0, 1, 1))
  recoveries <- allocate_vector("u32", 3L, 2L)   # 3 ticks x 2 nodes (time, node)

  sir_step(state, timer, nodeid, 4L, recoveries, 0L)

  expect_equal(state$values(), c(R, I, S, R))
  expect_equal(timer$values(), c(0L, 1L, 0L, 0L))
  m <- recoveries$values()                       # rows are ticks, columns are nodes
  expect_equal(m[1L, ], c(1, 1))                 # tick 0: one recovery in each node
  expect_equal(m[2L, ], c(0, 0))
  expect_equal(m[3L, ], c(0, 0))
})

test_that("sir_step writes each tick into its own recoveries column", {
  # Given the scenario above advanced one tick, where agent 2 (timer 2 -> 1) is
  #       still infectious
  # When a second sir_step runs for tick 1
  # Then agent 2 recovers into node 0's tick-1 column, while tick 0's column is
  #      left exactly as it was — confirming time-slices are independent.
  # Failure would mean sir_step writes the wrong column or clobbers earlier ticks.
  state  <- make("u8",  c(I, I, S, I))
  timer  <- make("u8",  c(1, 2, 0, 1))
  nodeid <- make("u16", c(0, 0, 1, 1))
  recoveries <- allocate_vector("u32", 3L, 2L)   # 3 ticks x 2 nodes

  sir_step(state, timer, nodeid, 4L, recoveries, 0L)
  sir_step(state, timer, nodeid, 4L, recoveries, 1L)

  expect_equal(state$values(), c(R, R, S, R))
  expect_equal(timer$values(), c(0L, 0L, 0L, 0L))
  m <- recoveries$values()
  expect_equal(m[1L, ], c(1, 1))                 # tick 0 unchanged
  expect_equal(m[2L, ], c(1, 0))                 # tick 1: agent 2 recovered in node 0
})

test_that("sir_step only processes the first `count` agents", {
  # Given 4 infectious agents but a count of 2
  # When sir_step runs
  # Then only agents 0 and 1 are considered: agent 0 (timer 1) recovers, agent 1
  #      (timer 5) stays I, and agents 2 and 3 are left infectious and untouched.
  # Failure would mean count is ignored and inactive slots are mutated.
  state  <- make("u8",  c(I, I, I, I))
  timer  <- make("u8",  c(1, 5, 1, 1))
  nodeid <- make("u16", c(0, 0, 1, 1))
  recoveries <- allocate_vector("u32", 1L, 2L)   # 1 tick x 2 nodes

  sir_step(state, timer, nodeid, 2L, recoveries, 0L)

  expect_equal(state$values(), c(R, I, I, I))    # agents 2,3 unchanged
  expect_equal(timer$values(), c(0L, 4L, 1L, 1L))
  expect_equal(recoveries$values(), c(1, 0))      # single tick -> plain vector
})

test_that("sir_step validates the tick index and the state/timer dtypes", {
  # Given out-of-range and wrong-dtype inputs
  # When sir_step is called
  # Then it errors rather than writing out of bounds or misreading the buffer.
  # Failure would risk silent memory corruption in the kernel.
  state  <- make("u8",  c(I, I))
  timer  <- make("u8",  c(1, 1))
  nodeid <- make("u16", c(0, 1))
  recoveries <- allocate_vector("u32", 2L, 2L)

  expect_error(sir_step(state, timer, nodeid, 2L, recoveries, 2L), "out of range")
  bad_state <- make("i32", c(I, I))             # state must be u8
  expect_error(sir_step(bad_state, timer, nodeid, 2L, recoveries, 0L), "u8")
})
