# Tests for bincount_where(): a predicate-filtered, count-aware bincount that answers
# flexible per-group agent queries (e.g. "exposed by node", "under-fives by node")
# directly on the Columns. Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]

test_that("bincount_where tallies agents in a state by node (ad-hoc vector mode)", {
  # Given 6 agents across 2 nodes with states S E E I E S
  # When we count those in state E by node
  # Then node 0 has 1 exposed and node 1 has 2 — matching a manual tabulate of the
  #      E agents' node ids. Failure would mean the predicate filter or grouping is wrong.
  state  <- allocate_scalar("u8",  6L); state$set(c(S, E, E, I, E, S))
  nodeid <- allocate_scalar("u16", 6L); nodeid$set(c(0, 0, 1, 1, 1, 0))

  got <- bincount_where(nodeid, 2L, state, "eq", E, count = 6L)

  expect_equal(got, c(1L, 2L))
  expect_equal(got, as.integer(tabulate(nodeid$values()[state$values() == E] + 1L, 2L)))
})

test_that("bincount_where finds under-fives by node via a dob > threshold predicate", {
  # Given agents with date-of-birth dob = -age (so younger = larger dob), under-five
  #       means tick - dob < 5*365, i.e. dob > tick - 5*365
  # When we count, per node, agents with dob greater than the threshold
  # Then only the genuinely young agents are counted per node. Failure would mean the
  #      comparison direction or the age arithmetic is wrong.
  tick <- 0L; thresh <- tick - 5L * 365L          # = -1825
  dob    <- allocate_scalar("i32", 5L); dob$set(c(-100, -2000, -50, -3000, -1800))  # ages 100,2000,50,3000,1800 days
  nodeid <- allocate_scalar("u16", 5L); nodeid$set(c(0, 0, 1, 1, 1))

  got <- bincount_where(nodeid, 2L, dob, "gt", thresh, count = 5L)

  # under five (dob > -1825): agents 1 (-100, node 0), 3 (-50, node 1), 5 (-1800, node 1)
  expect_equal(got, c(1L, 2L))
})

test_that("bincount_where scans only the live prefix `count`", {
  # Given a capacity-6 column with only the first 3 agents active
  # When we restrict the scan to count = 3
  # Then trailing (inactive) slots are ignored even if they would match. Failure would
  #      mean inactive/reserved slots leak into the tally.
  state  <- allocate_scalar("u8",  6L); state$set(c(E, E, S, E, E, E))   # slots 4-6 are E
  nodeid <- allocate_scalar("u16", 6L); nodeid$set(rep(0L, 6L))

  expect_equal(bincount_where(nodeid, 1L, state, "eq", E, count = 3L), 2L)  # only first 3 scanned
  expect_equal(bincount_where(nodeid, 1L, state, "eq", E, count = 6L), 5L)  # all 6 scanned
})

test_that("bincount_where writes into a report Column slice (model loop mode)", {
  # Given a 2-D (ticks x nodes) report Column
  # When bincount_where writes a tick's exposed-by-node into slice `slot`
  # Then that row holds the per-node counts and other rows are untouched — the pattern
  #      a model loop uses to record a per-tick, per-node series without reallocating.
  state  <- allocate_scalar("u8",  4L); state$set(c(E, E, S, E))
  nodeid <- allocate_scalar("u16", 4L); nodeid$set(c(0, 1, 1, 1))
  report <- allocate_vector("i32", 3L, 2L)         # 3 ticks x 2 nodes

  res <- bincount_where(nodeid, 2L, state, "eq", E, count = 4L, counts = report, slot = 1L)

  expect_null(res)
  expect_equal(report$values()[2L, ], c(1L, 2L))   # slot 1 -> row 2
  expect_true(all(report$values()[1L, ] == 0L))    # untouched
})

test_that("bincount_where rejects an unknown comparison op", {
  # Given an invalid op
  # When bincount_where is called
  # Then it errors rather than silently mis-counting.
  state  <- allocate_scalar("u8",  3L); state$set(c(E, E, S))
  nodeid <- allocate_scalar("u16", 3L); nodeid$set(rep(0L, 3L))
  expect_error(bincount_where(nodeid, 1L, state, "approx", E, count = 3L), "op")
})

test_that("bincount_where implements all six comparison ops correctly", {
  # Given a single group with properties 1..6 and a threshold of 3
  # When each of the six ops (eq/ne/lt/le/gt/ge) is applied
  # Then each count matches R's own comparison oracle. This pins the Rust op-parser against
  #      a swap (e.g. le/lt or ge/gt) — only "eq"/"gt" were previously exercised, yet
  #      run_model relies on "ne" and the API advertises all six. Failure means an op is
  #      mis-mapped and would silently mis-count.
  vals <- 1:6
  prop  <- allocate_scalar("i32", 6L); prop$set(vals)
  group <- allocate_scalar("u16", 6L); group$set(rep(0L, 6L))
  thr <- 3
  for (op in c("eq", "ne", "lt", "le", "gt", "ge")) {
    oracle <- switch(op,
      eq = sum(vals == thr), ne = sum(vals != thr),
      lt = sum(vals <  thr), le = sum(vals <= thr),
      gt = sum(vals >  thr), ge = sum(vals >= thr))
    expect_equal(bincount_where(group, 1L, prop, op, thr, count = 6L), as.integer(oracle),
                 info = paste("op", op))
  }
})

test_that("bincount_where requires an explicit count (no capacity default)", {
  # Given a call that omits count
  # When bincount_where / bincount_where_wt is invoked
  # Then it errors rather than defaulting to the Column's capacity (which would tally
  #      reserved/inactive slots in an over-allocated population). Failure would re-introduce
  #      the silent over-count.
  state  <- allocate_scalar("u8",  3L); state$set(c(I, I, S))
  nodeid <- allocate_scalar("u16", 3L); nodeid$set(rep(0L, 3L))
  wt     <- allocate_scalar("f64", 3L); wt$set(rep(1, 3L))
  expect_error(bincount_where(nodeid, 1L, state, "eq", I), "count")
  expect_error(bincount_where_wt(nodeid, 1L, state, "eq", I, wt), "count")
})

# ── bincount_where_wt: weighted + predicate-filtered ─────────────────────────────

test_that("bincount_where_wt sums a weight over matching agents by node (ad-hoc mode)", {
  # Given 5 agents (I I S I E) across 2 nodes with per-agent weights
  # When we sum weights of the infectious (state == I) by node
  # Then each node gets the sum of its infectious agents' weights only — matching a
  #      manual weighted tally. Failure would mean the predicate or weight sum is wrong.
  state  <- allocate_scalar("u8",  5L); state$set(c(I, I, S, I, E))
  nodeid <- allocate_scalar("u16", 5L); nodeid$set(c(0, 0, 1, 1, 1))
  shed   <- allocate_scalar("f64", 5L); shed$set(c(1.0, 0.5, 9, 2.0, 9))   # 9s on non-I

  got <- bincount_where_wt(nodeid, 2L, state, "eq", I, shed, count = 5L)

  expect_equal(got, c(1.5, 2.0))                      # node 0: 1.0+0.5; node 1: 2.0
  # equals the count-only result when every weight is 1:
  ones <- allocate_scalar("f64", 5L); ones$set(rep(1, 5L))
  expect_equal(bincount_where_wt(nodeid, 2L, state, "eq", I, ones, count = 5L),
               as.numeric(bincount_where(nodeid, 2L, state, "eq", I, count = 5L)))
})

test_that("bincount_where_wt writes into a report slice and honours the live prefix", {
  # Given a 2-D report and a weight column, with only the first 3 agents active
  # When bincount_where_wt writes node sums into slice `slot` over count = 3
  # Then trailing slots are ignored and the targeted slice holds the sums.
  state   <- allocate_scalar("u8",  5L); state$set(c(I, I, I, I, I))      # all I
  nodeid  <- allocate_scalar("u16", 5L); nodeid$set(c(0, 1, 1, 1, 1))
  weights <- allocate_scalar("f64", 5L); weights$set(c(2, 3, 4, 99, 99))  # slots 4-5 inactive
  report  <- allocate_vector("f64", 2L, 2L)

  res <- bincount_where_wt(nodeid, 2L, state, "eq", I, weights,
                           count = 3L, counts = report, slot = 1L)

  expect_null(res)
  expect_equal(report$values()[2L, ], c(2, 7))        # node 0: 2; node 1: 3+4 (slots 4-5 skipped)
  expect_true(all(report$values()[1L, ] == 0))
})
