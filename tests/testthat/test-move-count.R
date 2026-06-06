# Tests for move_count(): applies a per-node transition count returned by a kernel to the
# census at column tick+1 (from -= counts, to += counts), with either side NULL-skippable
# (one-sided decrement for deaths, increment for births). Written given-when-then.

test_that("move_count applies the delta to both compartments at tick+1", {
  # Given a from (10, 20) and to (1, 2) census, both carried to column 1
  # When move_count moves (3, 4) at tick 0
  # Then column 1 of from becomes (7, 16) and to becomes (4, 6); column 0 is untouched.
  from <- allocate_vector("i32", 2L, 2L); from$set(c(10L, 20L, 10L, 20L))
  to   <- allocate_vector("i32", 2L, 2L); to$set(c(1L, 2L, 1L, 2L))

  move_count(from, to, c(3L, 4L), 0L)

  expect_equal(from$col(1L), c(7L, 16L))
  expect_equal(to$col(1L),   c(4L, 6L))
  expect_equal(from$col(0L), c(10L, 20L))   # source tick untouched
})

test_that("move_count with a NULL side is one-sided (deaths / births)", {
  # Given a single census buffer carried to column 1
  # When move_count decrements only (to = NULL) and then increments only (from = NULL)
  # Then each one-sided update lands on column 1 as expected.
  col <- allocate_vector("i32", 2L, 2L); col$set(c(10L, 20L, 10L, 20L))

  move_count(col, NULL, c(2L, 3L), 0L)   # decrement only (e.g. a death)
  expect_equal(col$col(1L), c(8L, 17L))

  move_count(NULL, col, c(5L, 5L), 0L)   # increment only (e.g. a birth into M)
  expect_equal(col$col(1L), c(13L, 22L))
})
