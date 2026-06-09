# Tests for calc_capacity(): estimate the per-node agent capacity to preallocate for a
# population growing under a (time-varying) crude birth rate, with a safety-factor
# headroom. Ported from laser.core. Deterministic (no RNG). Written given-when-then.

test_that("zero birth rate yields a capacity equal to the initial population", {
  # Given no births anywhere over 100 steps and two nodes
  # When the capacity is estimated
  # Then growth is 1x and the safety factor adds nothing, so each capacity equals its
  #      initial population exactly. Failure would mean the growth/headroom math does
  #      not reduce to the identity at zero growth.
  br <- matrix(0, nrow = 100L, ncol = 2L)
  expect_equal(calc_capacity(br, c(1000, 250), safety_factor = 3), c(1000, 250))
})

test_that("a constant annual birth rate grows the population by ~ (1 + CBR/1000)", {
  # Given a constant CBR of 40 per 1,000 over exactly one year (365 steps) and no safety
  #       headroom (safety_factor = 0)
  # When the capacity is estimated for 1,000,000 initial
  # Then the estimate matches a full year of geometric growth, ~ initial*(1 + 40/1000)
  #      = 1,040,000, to within rounding/approximation. Failure would mean the daily-rate
  #      conversion or the exp(sum) growth is wrong.
  br <- matrix(40, nrow = 365L, ncol = 1L)
  cap <- calc_capacity(br, 1e6, safety_factor = 0)
  expect_equal(cap, 1.04e6, tolerance = 1e-5)   # ~1,040,002
})

test_that("the safety factor increases the estimate monotonically when growing", {
  # Given a positive birth rate
  # When the capacity is computed at safety factors 0, 1, 3
  # Then larger safety factors give larger capacities (more headroom), and all exceed
  #      the initial population. Failure would mean the safety multiplier is mis-signed
  #      or ignored.
  br <- matrix(30, nrow = 365L, ncol = 1L)
  c0 <- calc_capacity(br, 1e5, 0)
  c1 <- calc_capacity(br, 1e5, 1)
  c3 <- calc_capacity(br, 1e5, 3)
  expect_true(c0 < c1 && c1 < c3)
  expect_true(c0 > 1e5)
})

test_that("calc_capacity accepts a 2-D Column (values_map grid)", {
  # Given the SAME birth-rate grid expressed as a plain matrix and as a razer Column
  # When capacity is computed from each
  # Then the results are identical (the Column path just reads $values()). Failure would
  #      mean the Column branch reshapes or transposes the grid.
  nticks <- 50L; nnodes <- 3L
  br_mat <- matrix(rep(c(10, 20, 30), each = nticks), nrow = nticks, ncol = nnodes)
  br_col <- values_map(c(10, 20, 30), nticks, nnodes)   # per-node grid, same values
  pop <- c(1000, 2000, 3000)
  expect_equal(calc_capacity(br_col, pop), calc_capacity(br_mat, pop))
})

test_that("estimates above .Machine$integer.max are returned unclamped, with a warning", {
  # Given an enormous initial population and high sustained growth
  # When the capacity is estimated
  # Then the result is the true (large) double — NOT clamped — and a warning flags that
  #      it exceeds R's 32-bit integer max (the allocator's i32 count limit). Failure
  #      would mean silent truncation hiding an unallocatable capacity.
  br <- matrix(100, nrow = 3650L, ncol = 1L)        # max CBR for 10 years
  expect_warning(cap <- calc_capacity(br, 4e9, safety_factor = 6), "integer.max")
  expect_gt(cap, 2^32 - 1)                          # not clamped to uint32 max
  expect_lt(cap, 2^53)                              # still an exact whole-valued double
})

test_that("estimates within .Machine$integer.max produce no warning", {
  # Given a modest population and growth
  # When capacity is estimated
  # Then no warning is emitted (the result fits R's 32-bit integer range).
  # Failure would mean the warning threshold is mis-set and fires spuriously.
  br <- matrix(20, nrow = 365L, ncol = 2L)
  expect_no_warning(calc_capacity(br, c(1e6, 2e6)))
})

test_that("calc_capacity rejects invalid input", {
  # Given shape mismatches and out-of-range values
  # When calc_capacity is called
  # Then it errors rather than producing a nonsensical capacity.
  # Failure would mean bad inputs silently yield wrong preallocation.
  br <- matrix(20, nrow = 10L, ncol = 2L)
  expect_error(calc_capacity(br, c(100)), "must match")            # wrong node count
  expect_error(calc_capacity(br, c(-1, 100)), "non-negative")      # negative pop
  expect_error(calc_capacity(matrix(101, 5L, 1L), 100), "\\[0, 100\\]")  # CBR > 100
  expect_error(calc_capacity(matrix(-1, 5L, 1L), 100), "\\[0, 100\\]")   # CBR < 0
  expect_error(calc_capacity(br, c(1, 1), safety_factor = 7), "\\[0, 6\\]")  # bad SF
  expect_error(calc_capacity(c(1, 2, 3), 1), "2-D")                # not a 2-D grid
})

test_that("calc_capacity rejects NA / non-finite inputs with a clear message", {
  # Given NA or Inf in the birth-rate grid or the initial population
  # When calc_capacity is called
  # Then it errors with a clear 'finite' message rather than the cryptic base-R
  #      "missing value where TRUE/FALSE needed" that an unguarded comparison would raise.
  br <- matrix(30, 100L, 1L); br_na <- br; br_na[1L] <- NA
  expect_error(calc_capacity(br_na, 1e5), "finite")
  expect_error(calc_capacity(matrix(Inf, 5L, 1L), 100), "finite")
  expect_error(calc_capacity(br, NA), "finite")
})
