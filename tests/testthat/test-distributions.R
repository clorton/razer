# Distribution constructors (razer$normal, razer$constant) and their use as the
# infectious-period parameter of step_exposed_ei.
#
# The sampler RNG is thread-local (rand::thread_rng) and cannot be seeded from R,
# so stochastic tests assert distributional properties over large samples with
# wide margins rather than exact values.
#
# testthat idioms (for non-R readers): `test_that("desc", { ... })` is one test
# case. Matchers used here: `expect_true(cond)` / `expect_gt(a, b)` (a > b) /
# `expect_lt(a, b)` (a < b) / `expect_error(e)` (e raises — an extendr-surfaced
# Rust panic). `d <- dist_*(...)` constructs a sampler; `d$sample_n(k)` draws a
# length-k double vector.
# Vectorized R helpers used in assertions: `all(v)` is TRUE iff every element of a
# logical vector is TRUE; comparisons like `draws == 42`, `draws >= 3`,
# `draws > 0 & draws < 1` produce element-wise logical vectors; `mean`/`var`/`sd`/
# `median` are sample statistics; `abs(x) < tol` is the standard tolerance check.
# `&` is the element-wise (vectorized) AND.

# ── Constructors and $sample_one() ──────────────────────────────────────────────

test_that("dist_constant: every draw equals the supplied value", {
  # Given: a dist_constant(42) distribution
  # When:  sampled many times via $sample_one()
  # Then:  every draw is exactly 42
  # Failure would mean the degenerate distribution is not actually degenerate,
  # breaking its use as a fixed-duration drop-in.
  d <- dist_constant(42)
  draws <- d$sample_n(100L)
  # `draws == 42` is an element-wise comparison giving a logical vector; `all(...)`
  # is TRUE only if every element is TRUE.
  expect_true(all(draws == 42))
})

test_that("dist_normal: sample mean and variance match the parameters", {
  # Given: a dist_normal(mean = 10, variance = 4) distribution
  # When:  20 000 samples are drawn via $sample_one()
  # Then:  the sample mean is ~10 and the sample variance is ~4
  # Margins are wide enough to be robust to RNG variation (the standard error of
  # the mean for n = 20 000, sd = 2 is ~0.014, so [9.8, 10.2] is very safe).
  # Failure would indicate the mean/variance parameterization is wrong (e.g.
  # treating the second argument as a standard deviation).
  d <- dist_normal(10, 4)
  draws <- d$sample_n(20000L)
  # `expect_gt`/`expect_lt` are strict >/< assertions; bracketing the sample mean
  # in (9.8, 10.2) is the "approximately 10" check.
  expect_gt(mean(draws), 9.8)
  expect_lt(mean(draws), 10.2)
  expect_gt(var(draws), 3.0)
  expect_lt(var(draws), 5.0)
})

test_that("dist_normal: zero variance collapses to the mean", {
  # Given: a dist_normal(mean = 5, variance = 0) distribution
  # When:  sampled repeatedly
  # Then:  every draw is exactly the mean (5)
  # This documents the boundary case that a variance of 0 is accepted and yields
  # a deterministic value, matching dist_constant(5).
  d <- dist_normal(5, 0)
  draws <- d$sample_n(50L)
  expect_true(all(draws == 5))
})

test_that("dist_normal: negative variance is rejected", {
  # Given: a request for a normal distribution with negative variance
  # When:  dist_normal() is called
  # Then:  an error is raised (variance must be non-negative)
  # Failure would allow nonsensical parameters to construct an invalid sampler.
  expect_error(dist_normal(10, -1))
})

test_that("dist_uniform: draws fall within [low, high) with the correct mean", {
  # Given: a dist_uniform(3, 8) distribution
  # When:  20 000 samples are drawn via $sample_one()
  # Then:  every draw is in [3, 8) and the sample mean is ~5.5 (= (low + high) / 2)
  # Failure would indicate the bounds are not respected or are swapped.
  d <- dist_uniform(3, 8)
  draws <- d$sample_n(20000L)
  expect_true(all(draws >= 3))
  expect_true(all(draws < 8))
  expect_gt(mean(draws), 5.3)
  expect_lt(mean(draws), 5.7)
})

test_that("dist_uniform: a non-increasing interval is rejected", {
  # Given: degenerate or inverted bounds
  # When:  dist_uniform() is called
  # Then:  an error is raised (requires high > low)
  expect_error(dist_uniform(5, 5))
  expect_error(dist_uniform(8, 3))
})

test_that("dist_gamma: sample mean and variance match the shape-scale parameterization", {
  # Given: a dist_gamma(shape = 2, scale = 3) distribution (mean = 6, variance = 18)
  # When:  40 000 samples are drawn
  # Then:  all draws are strictly positive, the sample mean is ~6 and the sample
  #        variance is ~18
  # Margins are wide relative to the sampling error at n = 40 000. Failure would
  # indicate the (shape, scale) parameterization is wrong (e.g. shape-rate swap).
  d <- dist_gamma(2, 3)
  draws <- d$sample_n(40000L)
  expect_true(all(draws > 0))
  expect_gt(mean(draws), 5.6)
  expect_lt(mean(draws), 6.4)
  expect_gt(var(draws), 15)
  expect_lt(var(draws), 21)
})

test_that("dist_gamma: non-positive shape or scale is rejected", {
  # Given: invalid shape or scale parameters
  # When:  dist_gamma() is called
  # Then:  an error is raised (both must be positive)
  expect_error(dist_gamma(0, 3))
  expect_error(dist_gamma(2, -1))
})

test_that("dist_poisson: draws are non-negative integers with mean and variance ~ lambda", {
  # Given: a dist_poisson(5) distribution (mean = variance = 5)
  # When:  40 000 samples are drawn
  # Then:  all draws are non-negative integers, with sample mean ~5 and variance ~5
  # Failure would indicate the rate is misapplied or draws are not integer counts.
  d <- dist_poisson(5)
  draws <- d$sample_n(40000L)
  expect_true(all(draws >= 0))
  expect_true(all(draws == round(draws)))   # integer-valued counts (draws are doubles)
  expect_gt(mean(draws), 4.8)
  expect_lt(mean(draws), 5.2)
  expect_gt(var(draws), 4.4)
  expect_lt(var(draws), 5.6)
})

test_that("dist_poisson: non-positive lambda is rejected", {
  # Given: an invalid rate
  # When:  dist_poisson() is called
  # Then:  an error is raised (lambda must be positive)
  expect_error(dist_poisson(0))
  expect_error(dist_poisson(-2))
})

test_that("dist_beta: draws lie in (0, 1) with the correct mean and variance", {
  # Given: a dist_beta(2, 5) distribution (mean 2/7, variance 10/392 ≈ 0.0255)
  # When:  40 000 samples are drawn
  # Then:  every draw is strictly inside (0, 1), the sample mean is ~0.286 and the
  #        sample variance is ~0.0255
  # Failure would indicate the two shape parameters are swapped or mis-applied.
  d <- dist_beta(2, 5)
  draws <- d$sample_n(40000L)
  # `draws > 0 & draws < 1` is element-wise (vectorized) AND across the sample.
  expect_true(all(draws > 0 & draws < 1))
  # `abs(estimate - expected) < tol` is the two-sided tolerance check.
  expect_lt(abs(mean(draws) - 2 / 7),   0.01)
  expect_lt(abs(var(draws)  - 10 / 392), 0.005)
})

test_that("dist_beta: non-positive shapes are rejected", {
  # Given: invalid shape parameters
  # When:  dist_beta() is called
  # Then:  an error is raised (both shapes must be positive)
  expect_error(dist_beta(0, 5))
  expect_error(dist_beta(2, -1))
})

test_that("dist_exp: draws are positive with mean and sd equal to 1/rate", {
  # Given: a dist_exp(0.5) distribution (mean 2, variance 4, sd 2)
  # When:  40 000 samples are drawn
  # Then:  all draws are positive and the sample mean and sd are both ~2
  # Equal mean and sd is the exponential signature; both ~1/rate confirms the rate
  # is wired as a rate (not a mean/scale, which would give mean 0.5).
  d <- dist_exp(0.5)
  draws <- d$sample_n(40000L)
  expect_true(all(draws > 0))
  expect_lt(abs(mean(draws) - 2), 0.06)
  expect_lt(abs(sd(draws)   - 2), 0.06)
})

test_that("dist_exp: non-positive rate is rejected", {
  # Given: an invalid rate
  # When:  dist_exp() is called
  # Then:  an error is raised (rate must be positive)
  expect_error(dist_exp(0))
  expect_error(dist_exp(-1))
})

test_that("dist_logistic: symmetric about location with variance s^2 * pi^2 / 3", {
  # Given: a dist_logistic(4, 2) distribution (mean/median 4, variance 4*pi^2/3 ≈ 13.16)
  # When:  40 000 samples are drawn
  # Then:  the sample mean and median are ~4 and the sample variance is ~13.16
  # Failure would indicate the location/scale are mis-wired (e.g. scale ignored).
  d <- dist_logistic(4, 2)
  draws <- d$sample_n(40000L)
  expect_lt(abs(mean(draws)   - 4), 0.1)
  expect_lt(abs(median(draws) - 4), 0.1)
  expect_lt(abs(var(draws) - 4 * pi^2 / 3), 1.5)
})

test_that("dist_logistic: non-positive scale is rejected", {
  # Given: an invalid scale
  # When:  dist_logistic() is called
  # Then:  an error is raised (scale must be positive)
  expect_error(dist_logistic(4, 0))
  expect_error(dist_logistic(4, -2))
})

test_that("dist_lognormal: positive draws with median exp(meanlog) and mean exp(meanlog + sdlog^2/2)", {
  # Given: a dist_lognormal(0, 0.5) distribution (median exp(0) = 1, mean exp(0.125) ≈ 1.133)
  # When:  40 000 samples are drawn
  # Then:  all draws are positive, the sample median is ~1 and the sample mean is ~1.133
  # The median/mean split is the log-normal signature and confirms meanlog and sdlog
  # are treated as log-space parameters (matching R's qlnorm).
  d <- dist_lognormal(0, 0.5)
  draws <- d$sample_n(40000L)
  expect_true(all(draws > 0))
  expect_lt(abs(median(draws) - 1),               0.05)
  expect_lt(abs(mean(draws)   - exp(0.5^2 / 2)),   0.05)
})

test_that("dist_lognormal: negative sdlog is rejected", {
  # Given: an invalid (negative) log-space standard deviation
  # When:  dist_lognormal() is called
  # Then:  an error is raised (sdlog must be non-negative)
  expect_error(dist_lognormal(0, -0.5))
})

# ── Use as the infectious-period distribution of step_exposed_ei ─────────────────

# A local test fixture/factory. `function(args)` defines an anonymous function
# bound to `make_exposed`; `timer = 1L` is a default argument.
make_exposed <- function(n, timer = 1L) {
  # n exposed (state E) agents in a single node, each with the given timer so
  # that one step of step_exposed_ei expires the timer and forces E→I.
  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", 1L)   # E
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", timer)
  ppl
}

test_that("step_exposed_ei: constant distribution sets a fixed infectious timer", {
  # Given: 1 000 exposed agents whose timers expire this tick
  # When:  step_exposed_ei is run with a dist_constant(9) infectious-period distribution
  # Then:  all agents become infectious (state I) with timer exactly 9
  # Failure would mean the constant distribution is not being sampled as a fixed
  # duration on the E→I transition.
  ppl <- make_exposed(1000L, timer = 1L)
  # `inf_dist =` passes the distribution by keyword argument. The frame `ppl` is a
  # reference (external pointer), so the step mutates it in place.
  step_exposed_ei(ppl, inf_dist = dist_constant(9))
  expect_true(all(ppl$state == 2L))
  expect_true(all(ppl$timer == 9L))
})

test_that("step_exposed_ei: normal distribution produces varied integer timers", {
  # Given: 20 000 exposed agents whose timers expire this tick
  # When:  step_exposed_ei is run with a dist_normal(mean = 12, variance = 9) period
  # Then:  all transition to I; their infectious timers are integers >= 1, vary
  #        across agents (sd > 0), and average near the requested mean of 12
  # This is the headline behavior: a per-agent stochastic infectious period drawn
  # from the supplied distribution. Failure indicates the distribution is not
  # actually sampled per agent (e.g. a single shared draw) or the step's
  # round-and-clamp to whole ticks is wrong.
  ppl <- make_exposed(20000L, timer = 1L)
  step_exposed_ei(ppl, inf_dist = dist_normal(12, 9))

  timers <- ppl$timer
  expect_true(all(ppl$state == 2L))
  # `as.integer` truncates to whole numbers; equality with the originals asserts
  # the durations are already whole ticks.
  expect_true(all(timers == as.integer(timers)))  # whole-tick durations
  expect_true(all(timers >= 1L))                   # clamped to a positive period
  expect_gt(sd(timers), 1)                         # genuinely stochastic spread
  expect_gt(mean(timers), 11)
  expect_lt(mean(timers), 13)
})

test_that("step_transmission_si: newly infected agents draw their timer from any family", {
  # Given: a single node where every susceptible is infected this tick (beta huge),
  #        seeded with infectious agents to drive the force of infection, and a
  #        dist_uniform(10, 20) infectious-period distribution
  # When:  one transmission step is run
  # Then:  every newly infected agent gets an integer timer within the support of
  #        the distribution (rounded dist_uniform(10, 20) lies in 10:20), and those
  #        timers vary across agents
  # This confirms the new distribution families flow through a transmission kernel
  # (not just the timer-decrement kernels) and are sampled per agent.
  n <- 10000L
  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", 0L)              # S
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)
  # Build the initial state vector in R, then write it to the frame. `;` separates
  # statements on one line; `seq_len(100L)` is 1:100 (but safe when the count is 0);
  # `sv[1:100] <- 2L` assigns into a slice; `2L` is the I (infectious) state code.
  sv <- rep(0L, n); sv[seq_len(100L)] <- 2L; ppl$state <- sv   # 100 seeded I

  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  state_before <- ppl$state
  step_transmission_si(ppl, nd, beta = 100.0, inf_dist = dist_uniform(10, 20), network = matrix(0, 1, 1))

  # `which(logical_vec)` returns the 1-based indices where the condition is TRUE
  # (here: was S, now I). `ppl$timer[newly_infected]` then gathers those agents' timers.
  newly_infected <- which(state_before == 0L & ppl$state == 2L)
  timers <- ppl$timer[newly_infected]
  expect_gt(length(newly_infected), 0L)
  expect_true(all(timers == as.integer(timers)))   # whole-tick durations
  expect_true(all(timers >= 10L & timers <= 20L))   # within rounded uniform support
  expect_gt(sd(timers), 0)                          # per-agent draws, not a shared value
})
