# Tests for AliasedDistribution: a Vose alias-method sampler over 0-based bin indices,
# weighted by integer counts (e.g. a population pyramid's per-age-band counts). Draws
# are 0-based. The RNG is thread-local (not R-seedable), so the sampling tests are
# statistical with generous tolerances at large sample sizes. Written given-when-then.

test_that("aliased_distribution samples bins in proportion to their counts", {
  # Given counts 10, 30, 60 (proportions 0.1, 0.3, 0.6)
  # When 200,000 bins are drawn
  # Then the empirical bin frequencies match the proportions to within 0.01, and the
  #      returned indices are 0-based (in 0..2). Failure would mean the alias table or
  #      the threshold comparison is wrong, biasing the generated population.
  d <- aliased_distribution(c(10, 30, 60))
  s <- d$sample_n(200000L)

  expect_equal(range(s), c(0L, 2L))                      # 0-based bin indices
  freqs <- tabulate(s + 1L, nbins = 3L) / length(s)      # +1L: 0-based -> 1-based
  expect_equal(freqs, c(0.1, 0.3, 0.6), tolerance = 0.01)
})

test_that("aliased_distribution with equal counts is uniform", {
  # Given four equal counts
  # When 100,000 bins are drawn
  # Then each bin is drawn ~25% of the time (no aliasing needed; every column is
  #      exactly full). Failure would mean the "exactly full" bins are mishandled.
  d <- aliased_distribution(c(5, 5, 5, 5))
  freqs <- tabulate(d$sample_n(100000L) + 1L, nbins = 4L) / 100000
  expect_equal(freqs, rep(0.25, 4L), tolerance = 0.01)
})

test_that("a single-bin distribution always returns bin 0", {
  # Given one bin
  # When samples are drawn
  # Then every draw is bin 0. Failure would mean a degenerate table mis-indexes.
  d <- aliased_distribution(c(42))
  expect_equal(d$n_bins(), 1L)
  expect_true(all(d$sample_n(1000L) == 0L))
  expect_equal(d$sample_one(), 0L)
})

test_that("aliased_distribution exposes its bin count and total weight", {
  # Given counts whose sum exceeds R's integer range
  # When the accessors are queried
  # Then n_bins is the bin count and total is the (double) sum. Failure would indicate
  #      an overflow in the integer-scaled alias construction.
  d <- aliased_distribution(c(2e9, 2e9))                 # sum 4e9 > .Machine$integer.max
  expect_equal(d$n_bins(), 2L)
  expect_equal(d$total(), 4e9)
  expect_equal(length(d$alias()), 2L)
  expect_equal(length(d$probs()), 2L)
})

test_that("aliased_distribution rejects invalid counts", {
  # Given degenerate or invalid count vectors
  # When the constructor is called
  # Then it errors rather than building an unusable sampler.
  # Failure would risk a divide-by-zero or a biased/looping construction.
  expect_error(aliased_distribution(numeric(0)), "at least one bin")
  expect_error(aliased_distribution(c(0, 0, 0)), "positive total")
  expect_error(aliased_distribution(c(1, -2, 3)), "non-negative")
})
