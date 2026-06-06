# Tests for import_infections(): activates RESERVED agent slots (those past the live
# `count`, in capacity-sized property arrays) as new infectious cases drawn from a
# schedule (parallel tick / node / count vectors). For the entries matching the current
# tick it sets each new agent's state to Infectious, timer from a duration Distribution,
# and nodeid to its node; adds the per-node imports to the I census at tick+1 and to the
# importations flow at tick; and returns the grown live count. Written given-when-then.

states <- laser_states()
S <- states[["S"]]; I <- states[["I"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }

test_that("import_infections activates reserved slots and updates census + flow", {
  # Given 4 live agents in capacity-6 arrays (2 reserved slots), an I census and an
  #       importations flow over 3 ticks / 2 nodes, and a schedule importing 1 case
  #       into node 0 at tick 0 and 1 into node 1 at tick 1, with a constant timer 7
  # When import_infections runs for tick 0
  # Then slot 4 (the first reserved) becomes Infectious in node 0 with timer 7; the
  #      live count grows to 5; I census at tick+1 gains (1, 0); the importations flow
  #      at tick 0 is (1, 0); and tick 1's entry is left alone (processed on its tick).
  # Failure would mean imports land in the wrong slot/node/tick or miscount the census.
  state  <- mk("u8",  c(I, I, I, I, S, S))     # slots 4,5 reserved (state S, inactive)
  timer  <- mk("u16",  c(3, 3, 3, 3, 0, 0))
  nodeid <- mk("u16", c(0, 0, 1, 1, 0, 0))
  i_count <- allocate_vector("i32", 3L, 2L)
  importations <- allocate_vector("i32", 2L, 2L)
  duration <- dist_constant(7)

  new_count <- import_infections(
    state, timer, nodeid, 4L, i_count, importations,
    sched_tick  = c(0L, 1L),
    sched_node  = c(0L, 1L),
    sched_count = c(1L, 1L),
    duration, 0L)

  expect_equal(new_count, 5L)                       # one agent imported at tick 0
  sv <- state$values()
  expect_equal(sv[5L], I)                           # slot 4 (0-based) now infectious
  expect_equal(timer$values()[5L], 7L)              # timer from the duration
  expect_equal(nodeid$values()[5L], 0L)             # placed in node 0
  expect_equal(sv[6L], S)                           # slot 5 still reserved/inactive
  expect_equal(i_count$values()[2L, ], c(1L, 0L))   # I census at tick+1
  expect_equal(importations$values()[1L, ], c(1L, 0L))  # flow at tick 0
})

test_that("import_infections imports nothing on a tick with no schedule entries", {
  # Given a schedule whose only entry is at tick 1
  # When import_infections runs for tick 0
  # Then the count is unchanged, no slot is touched, and the census/flow stay zero.
  # Failure would mean it imports on the wrong tick.
  state  <- mk("u8",  c(I, S)); timer <- mk("u16", c(2, 0)); nodeid <- mk("u16", c(0, 0))
  i_count <- allocate_vector("i32", 3L, 1L)
  importations <- allocate_vector("i32", 2L, 1L)

  new_count <- import_infections(
    state, timer, nodeid, 1L, i_count, importations,
    c(1L), c(0L), c(1L), dist_constant(5), 0L)

  expect_equal(new_count, 1L)
  expect_equal(state$values(), c(I, S))
  expect_equal(i_count$values()[2L, ], 0)
  expect_equal(importations$values()[1L, ], 0)
})

test_that("import_infections aggregates multiple entries for the same tick", {
  # Given two schedule entries at tick 0 (2 into node 0, 3 into node 1) with 5 reserved
  #       slots
  # When import_infections runs for tick 0
  # Then 5 agents activate (count 0 -> 5), the I census at tick+1 is (2, 3), and the
  #      flow is (2, 3). Failure would mean per-tick entries are not summed per node.
  state  <- mk("u8",  rep(S, 5L)); timer <- mk("u16", rep(0, 5L)); nodeid <- mk("u16", rep(0, 5L))
  i_count <- allocate_vector("i32", 2L, 2L)
  importations <- allocate_vector("i32", 1L, 2L)

  new_count <- import_infections(
    state, timer, nodeid, 0L, i_count, importations,
    c(0L, 0L), c(0L, 1L), c(2L, 3L), dist_constant(4), 0L)

  expect_equal(new_count, 5L)
  expect_equal(sum(state$values() == I), 5L)
  expect_equal(i_count$values()[2L, ], c(2L, 3L))
  expect_equal(importations$values(), c(2L, 3L))    # single-tick buffer: plain vector
})

test_that("import_infections errors when imports would exceed capacity", {
  # Given 2 live agents in capacity-2 arrays (no reserved room) and a schedule asking
  #       to import 1 more
  # When import_infections runs
  # Then it errors rather than writing past the allocated buffer.
  # Failure would risk an out-of-bounds slot write.
  state  <- mk("u8",  c(I, I)); timer <- mk("u16", c(1, 1)); nodeid <- mk("u16", c(0, 0))
  i_count <- allocate_vector("i32", 2L, 1L)
  importations <- allocate_vector("i32", 1L, 1L)

  expect_error(
    import_infections(state, timer, nodeid, 2L, i_count, importations,
                      c(0L), c(0L), c(1L), dist_constant(3), 0L),
    "exceed capacity")
})
