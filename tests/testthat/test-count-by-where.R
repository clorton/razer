# Tests for count_by_where(): a predicate-filtered, count-aware bincount that answers
# flexible per-group agent queries (e.g. "exposed by node", "under-fives by node")
# directly on the Columns. Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]

test_that("count_by_where tallies agents in a state by node (ad-hoc vector mode)", {
  # Given 6 agents across 2 nodes with states S E E I E S
  # When we count those in state E by node
  # Then node 0 has 1 exposed and node 1 has 2 — matching a manual tabulate of the
  #      E agents' node ids. Failure would mean the predicate filter or grouping is wrong.
  state  <- allocate_scalar("u8",  6L); state$set(c(S, E, E, I, E, S))
  nodeid <- allocate_scalar("u16", 6L); nodeid$set(c(0, 0, 1, 1, 1, 0))

  got <- count_by_where(nodeid, 2L, state, "eq", E)

  expect_equal(got, c(1L, 2L))
  expect_equal(got, as.integer(tabulate(nodeid$values()[state$values() == E] + 1L, 2L)))
})

test_that("count_by_where finds under-fives by node via a dob > threshold predicate", {
  # Given agents with date-of-birth dob = -age (so younger = larger dob), under-five
  #       means tick - dob < 5*365, i.e. dob > tick - 5*365
  # When we count, per node, agents with dob greater than the threshold
  # Then only the genuinely young agents are counted per node. Failure would mean the
  #      comparison direction or the age arithmetic is wrong.
  tick <- 0L; thresh <- tick - 5L * 365L          # = -1825
  dob    <- allocate_scalar("i32", 5L); dob$set(c(-100, -2000, -50, -3000, -1800))  # ages 100,2000,50,3000,1800 days
  nodeid <- allocate_scalar("u16", 5L); nodeid$set(c(0, 0, 1, 1, 1))

  got <- count_by_where(nodeid, 2L, dob, "gt", thresh)

  # under five (dob > -1825): agents 1 (-100, node 0), 3 (-50, node 1), 5 (-1800, node 1)
  expect_equal(got, c(1L, 2L))
})

test_that("count_by_where scans only the live prefix `count`", {
  # Given a capacity-6 column with only the first 3 agents active
  # When we restrict the scan to count = 3
  # Then trailing (inactive) slots are ignored even if they would match. Failure would
  #      mean inactive/reserved slots leak into the tally.
  state  <- allocate_scalar("u8",  6L); state$set(c(E, E, S, E, E, E))   # slots 4-6 are E
  nodeid <- allocate_scalar("u16", 6L); nodeid$set(rep(0L, 6L))

  expect_equal(count_by_where(nodeid, 1L, state, "eq", E, count = 3L), 2L)  # only first 3 scanned
  expect_equal(count_by_where(nodeid, 1L, state, "eq", E),             5L)  # default: all 6
})

test_that("count_by_where writes into a report Column slice (model loop mode)", {
  # Given a 2-D (ticks x nodes) report Column
  # When count_by_where writes a tick's exposed-by-node into slice `slot`
  # Then that row holds the per-node counts and other rows are untouched — the pattern
  #      a model loop uses to record a per-tick, per-node series without reallocating.
  state  <- allocate_scalar("u8",  4L); state$set(c(E, E, S, E))
  nodeid <- allocate_scalar("u16", 4L); nodeid$set(c(0, 1, 1, 1))
  report <- allocate_vector("i32", 3L, 2L)         # 3 ticks x 2 nodes

  res <- count_by_where(nodeid, 2L, state, "eq", E, counts = report, slot = 1L)

  expect_null(res)
  expect_equal(report$values()[2L, ], c(1L, 2L))   # slot 1 -> row 2
  expect_true(all(report$values()[1L, ] == 0L))    # untouched
})

test_that("count_by_where rejects an unknown comparison op", {
  # Given an invalid op
  # When count_by_where is called
  # Then it errors rather than silently mis-counting.
  state  <- allocate_scalar("u8",  3L); state$set(c(E, E, S))
  nodeid <- allocate_scalar("u16", 3L); nodeid$set(rep(0L, 3L))
  expect_error(count_by_where(nodeid, 1L, state, "approx", E), "op")
})
