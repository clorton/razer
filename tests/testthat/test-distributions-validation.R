# Statistical validation of the Rust distribution samplers against R's own
# reference implementations (stats::qnorm / qunif / qgamma / qpois / dpois).
#
# Purpose: confirm that the parameters passed from R are wired into the Rust crate
# correctly — i.e. that `dist_normal(mean, variance)` really produces N(mean, sd =
# sqrt(variance)), that `dist_gamma(shape, scale)` uses the shape-SCALE (not
# shape-rate) parameterization, that `dist_poisson(lambda)` has mean = variance =
# lambda, etc. Each test draws a large empirical sample in a single batch
# (`$sample_n()`) and compares moments and quantiles to the known theoretical
# values. A wiring mistake (argument swap, sd-vs-variance, rate-vs-scale) shifts
# the empirical distribution far enough to fail with overwhelming margin.
#
# The RNG (rand::thread_rng) is not seedable from R, so tolerances are set wide
# relative to the Monte-Carlo standard error at N = 1e6 (many standard errors), so
# the tests are robust to ordinary sampling variation but still tight enough to
# catch a misparameterization.

N <- 1000000L

# Compare empirical quantiles to a theoretical quantile function at several probs.
expect_quantiles_close <- function(draws, qfun, probs, tol) {
  emp  <- as.numeric(stats::quantile(draws, probs, names = FALSE, type = 7))
  theo <- qfun(probs)
  expect_lt(max(abs(emp - theo)), tol)
}

PROBS <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)

# ── sampler plumbing ─────────────────────────────────────────────────────────────

test_that("sample_n / sample_one return doubles of the expected length", {
  # Given: any distribution
  # When:  sample_one() and sample_n(k) are called
  # Then:  sample_one() yields a single finite double; sample_n(k) yields a length-k
  #        double vector; sample_n(0) yields a length-0 vector
  # Failure would indicate the batch sampler is mis-sized or mistyped, invalidating
  # every downstream validation below.
  d <- dist_normal(0, 1)
  one <- d$sample_one()
  expect_type(one, "double")
  expect_length(one, 1L)
  expect_true(is.finite(one))

  expect_length(d$sample_n(1000L), 1000L)
  expect_length(d$sample_n(0L), 0L)
})

# ── normal: dist_normal(mean, variance) == N(mean, sd = sqrt(variance)) ───────────

test_that("dist_normal is wired as N(mean, variance) and matches qnorm", {
  # Given: dist_normal(10, 4)  ->  mean 10, variance 4, sd 2
  # When:  one million draws are taken
  # Then:  sample mean ~10, sample sd ~2 (NOT 4 — the arg is variance, not sd), and
  #        empirical quantiles match qnorm(p, 10, 2)
  # The sd check is the key wiring assertion: if variance were treated as sd, the
  # sample sd would be ~4 and both it and the quantiles would fail by a wide margin.
  draws <- dist_normal(10, 4)$sample_n(N)
  expect_lt(abs(mean(draws) - 10),  0.05)
  expect_lt(abs(sd(draws)   -  2),  0.02)
  expect_quantiles_close(draws, function(p) qnorm(p, mean = 10, sd = 2), PROBS, tol = 0.05)
})

# ── uniform: dist_uniform(low, high) == U[low, high) ──────────────────────────────

test_that("dist_uniform is wired as U(low, high) and matches qunif", {
  # Given: dist_uniform(2, 10)  ->  support [2, 10), mean 6, variance (8^2)/12 = 5.333
  # When:  one million draws are taken
  # Then:  all draws lie in [2, 10), sample mean ~6, sample variance ~5.333, and
  #        empirical quantiles match qunif(p, 2, 10)
  draws <- dist_uniform(2, 10)$sample_n(N)
  expect_gte(min(draws), 2)
  expect_lt(max(draws), 10)
  expect_lt(abs(mean(draws) - 6),         0.02)
  expect_lt(abs(var(draws)  - 64 / 12),   0.05)
  expect_quantiles_close(draws, function(p) qunif(p, min = 2, max = 10), PROBS, tol = 0.05)
})

# ── gamma: dist_gamma(shape, scale) uses the shape-SCALE parameterization ──────────

test_that("dist_gamma is wired with shape-scale parameterization and matches qgamma", {
  # Given: dist_gamma(3, 2)  ->  mean shape*scale = 6, variance shape*scale^2 = 12
  # When:  one million draws are taken
  # Then:  all draws > 0, sample mean ~6, sample variance ~12, and empirical
  #        quantiles match qgamma(p, shape = 3, scale = 2)
  # The mean check distinguishes scale from rate: a shape-RATE reading of (3, 2)
  # would give mean 3/2 = 1.5, failing immediately.
  draws <- dist_gamma(3, 2)$sample_n(N)
  expect_gt(min(draws), 0)
  expect_lt(abs(mean(draws) - 6),   0.05)
  expect_lt(abs(var(draws)  - 12),  0.3)
  expect_quantiles_close(draws, function(p) qgamma(p, shape = 3, scale = 2), PROBS, tol = 0.2)
})

# ── poisson: dist_poisson(lambda) has mean = var = lambda; PMF matches dpois ───────

test_that("dist_poisson is wired with rate lambda and matches dpois", {
  # Given: dist_poisson(7)  ->  mean 7, variance 7, non-negative integer counts
  # When:  one million draws are taken
  # Then:  all draws are non-negative integers, sample mean ~7 and variance ~7, and
  #        the empirical PMF over counts 0..18 matches dpois(0:18, 7)
  # Equality of mean and variance (both ~lambda) is the Poisson signature; matching
  # the full low-order PMF confirms the rate parameter is wired correctly.
  draws <- dist_poisson(7)$sample_n(N)
  expect_true(all(draws >= 0))
  expect_true(all(draws == round(draws)))
  expect_lt(abs(mean(draws) - 7),  0.05)
  expect_lt(abs(var(draws)  - 7),  0.1)

  ks  <- 0:18
  emp <- vapply(ks, function(k) mean(draws == k), numeric(1))
  expect_lt(max(abs(emp - dpois(ks, lambda = 7))), 0.005)
})
