# Tests for calc_capacity_cdr(): the mortality-aware capacity estimate for populations
# whose dead slots are reclaimed periodically with squash(). It bounds the PEAK LIVING
# population (net births minus a conservatively underestimated death rate), not the
# cumulative number ever born. Deterministic, so exact assertions. Written given-when-then.

# Closed-form reference for the documented formula.
ref_cap <- function(cbr, cdr, nsteps, pop, sf) {
  lb <- (1 + cbr / 1000)^(1 / 365) - 1
  ld <- (1 + cdr / 1000)^(1 / 365) - 1
  round(pop * exp(nsteps * lb - (1 / (1 + sf)) * nsteps * ld))
}

test_that("calc_capacity_cdr matches its closed-form net-growth formula", {
  # Given 1,000,000 people over 100 years at CBR 30 / CDR 15, safety_factor 1
  # When calc_capacity_cdr estimates the capacity
  # Then it equals pop * exp(sum lambda_b - 0.5 * sum lambda_d) — the peak-living bound with
  #      deaths credited at 1/(1+sf). Failure would mean the formula/rate conversion is wrong.
  nsteps <- 100L * 365L
  br <- matrix(30, nrow = nsteps, ncol = 1)
  dr <- matrix(15, nrow = nsteps, ncol = 1)

  got <- calc_capacity_cdr(br, dr, initial_pop = 1e6, safety_factor = 1)

  expect_equal(got, ref_cap(30, 15, nsteps, 1e6, 1))
})

test_that("calc_capacity_cdr needs far fewer slots than the cumulative-births bound", {
  # Given the same century-long CBR 30 / CDR 15 scenario
  # When compared to calc_capacity() (which assumes NO reclaim, so counts all births)
  # Then the squash-aware estimate is much smaller yet still above the bare net-growth peak
  #      (sf = 0) — i.e. reclaiming dead slots genuinely saves memory. Failure would mean
  #      the mortality credit isn't reducing the bound.
  nsteps <- 100L * 365L
  br <- matrix(30, nrow = nsteps, ncol = 1)
  dr <- matrix(15, nrow = nsteps, ncol = 1)

  cdr_cap   <- calc_capacity_cdr(br, dr, 1e6, safety_factor = 1)
  births_cap <- calc_capacity(br, 1e6, safety_factor = 1)
  net_peak  <- calc_capacity_cdr(br, dr, 1e6, safety_factor = 0)   # full death credit = bare net

  expect_lt(cdr_cap, births_cap)        # squash saves slots vs the no-reclaim bound
  expect_gt(cdr_cap, net_peak)          # but reserves headroom above the bare net peak
  expect_equal(net_peak, ref_cap(30, 15, nsteps, 1e6, 0))
})

test_that("calc_capacity_cdr increases with the safety factor (more death headroom)", {
  # Given increasing safety factors
  # When the death credit shrinks (1/(1+sf))
  # Then the estimated capacity grows monotonically. Failure would mean the safety knob is
  #      inverted or ignored.
  nsteps <- 50L * 365L
  br <- matrix(25, nrow = nsteps, ncol = 1)
  dr <- matrix(10, nrow = nsteps, ncol = 1)
  caps <- vapply(c(0, 0.5, 1, 2, 4), function(sf) calc_capacity_cdr(br, dr, 5e5, sf), numeric(1))
  expect_true(all(diff(caps) > 0))
})

test_that("calc_capacity_cdr handles multiple nodes and Column inputs", {
  # Given per-node rate grids supplied both as matrices and as values_map Columns
  # When estimating capacity for two nodes
  # Then each node is estimated independently and the Column path agrees with the matrix
  #      path. Failure would mean node handling or the Column coercion is wrong.
  nsteps <- 10L * 365L; nn <- 2L
  br <- matrix(rep(c(30, 20), each = nsteps), ncol = nn)
  dr <- matrix(rep(c(15, 18), each = nsteps), ncol = nn)
  m_caps <- calc_capacity_cdr(br, dr, c(1e5, 2e5), safety_factor = 1)
  expect_length(m_caps, nn)
  expect_equal(m_caps[1L], ref_cap(30, 15, nsteps, 1e5, 1))
  expect_equal(m_caps[2L], ref_cap(20, 18, nsteps, 2e5, 1))

  c_caps <- calc_capacity_cdr(values_map(c(30, 20), nsteps, nn),
                              values_map(c(15, 18), nsteps, nn), c(1e5, 2e5), 1)
  expect_equal(c_caps, m_caps)
})

test_that("calc_capacity_cdr validates its inputs", {
  # Given contract violations
  # When calc_capacity_cdr is called
  # Then each raises an informative error rather than returning nonsense.
  br <- matrix(30, 100, 2); dr <- matrix(15, 100, 2)
  expect_error(calc_capacity_cdr(br, matrix(15, 100, 3), c(1, 1), 1), "same shape")
  expect_error(calc_capacity_cdr(br, dr, c(1, 1, 1), 1), "must match")
  expect_error(calc_capacity_cdr(br, dr, c(-1, 1), 1), "non-negative")
  expect_error(calc_capacity_cdr(matrix(200, 100, 2), dr, c(1, 1), 1), "0, 100")
  expect_error(calc_capacity_cdr(br, dr, c(1, 1), safety_factor = 7), "0, 6")
})
