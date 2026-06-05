# Tests for the migration-network models ported from laser-core (gravity,
# radiation, stouffer, competing_destinations) and the row_normalizer helper.
# Written given-when-then. These models are deterministic, so the expected values
# below are regression anchors computed from laser-core's reference numpy
# implementation on the same 4-node scenario.

# A small, hand-traceable 4-node scenario reused across tests: populations and a
# symmetric distance matrix (`byrow = TRUE` fills row-by-row, the human reading
# order; R stores it column-major internally).
pops <- c(1000, 500, 200, 800)
D <- matrix(c(0, 2, 3, 5,
              2, 0, 1, 4,
              3, 1, 0, 2,
              5, 4, 2, 0), nrow = 4L, byrow = TRUE)

test_that("gravity matches the laser-core reference and is symmetric when a == b", {
  # Given the 4-node scenario and gravity parameters k=0.5, a=b=1, c=2
  # When the gravity network is computed
  # Then it matches laser-core's reference values, has a zero diagonal, and is
  #      symmetric (because a == b makes p_i^a p_j^b symmetric in i, j).
  # Failure would mean the gravity formula or the column-major layout is wrong.
  g <- gravity(pops, D, k = 0.5, a = 1, b = 1, c = 2)

  expect_equal(diag(g), c(0, 0, 0, 0))
  expect_equal(g, t(g))                       # a == b  =>  symmetric
  expect_equal(g[1L, 2L], 0.5 * 1000 * 500 / 2^2)   # 62500
  expect_equal(g[1L, 3L], 0.5 * 1000 * 200 / 3^2)   # ~11111.11
})

test_that("radiation matches the laser-core reference (include_home = FALSE)", {
  # Given the 4-node scenario, k=1, include_home=FALSE
  # When the radiation network is computed
  # Then every entry matches laser-core's reference radiation output.
  # Failure would mean the distance ranking, the as-close-or-closer cumulative
  #      sum, or the home-exclusion is wrong.
  expected <- matrix(c(
    0.000000, 0.166667, 0.061920, 0.096970,
    0.108932, 0.000000, 0.158730, 0.048485,
    0.022857, 0.119048, 0.000000, 0.046377,
    0.091429, 0.133333, 0.133333, 0.000000),
    nrow = 4L, byrow = TRUE)

  r <- radiation(pops, D, k = 1, include_home = FALSE)

  expect_equal(diag(r), c(0, 0, 0, 0))
  expect_equal(r, expected, tolerance = 1e-5)
})

test_that("radiation include_home = TRUE yields smaller weights than FALSE", {
  # Given the same scenario under both include_home settings
  # When both radiation networks are computed
  # Then including the home population in the cumulative denominator strictly
  #      lowers every off-diagonal weight (larger denominator).
  # Failure would mean the include_home flag is ignored or inverted.
  r_false <- radiation(pops, D, k = 1, include_home = FALSE)
  r_true  <- radiation(pops, D, k = 1, include_home = TRUE)

  off <- D > 0                                # off-diagonal mask
  expect_true(all(r_true[off] < r_false[off]))
})

test_that("stouffer matches the laser-core reference", {
  # Given the 4-node scenario, k=1, a=b=1, include_home=FALSE
  # When the stouffer network is computed
  # Then it matches laser-core's reference output and has a zero diagonal.
  # Failure would mean the nearest-node (home) skip or the ratio^b term is wrong.
  expected <- matrix(c(
      0.000000, 1000.0000, 285.71429, 533.33333,
    416.666667,    0.0000, 500.00000, 200.00000,
     86.956522,  200.0000,   0.00000, 123.07692,
    470.588235,  571.4286, 800.00000,   0.00000),
    nrow = 4L, byrow = TRUE)

  s <- stouffer(pops, D, k = 1, a = 1, b = 1, include_home = FALSE)

  expect_equal(diag(s), c(0, 0, 0, 0))
  expect_equal(s, expected, tolerance = 1e-4)
})

test_that("competing_destinations matches the laser-core reference", {
  # Given the 4-node scenario, gravity params k=0.5, a=b=1, c=2 and delta=0.5
  # When the competing-destinations network is computed
  # Then it matches laser-core's reference output.
  # Failure would mean the accessibility (competition) term or its k!=i, k!=j
  #      exclusions are wrong.
  expected <- matrix(c(
         0.0, 988211.8, 293972.4, 144222.1,
    460223.4,      0.0, 881917.1, 118585.4,
    139221.8, 866025.4,      0.0, 168819.4,
    194136.3, 265165.0, 494413.2,      0.0),
    nrow = 4L, byrow = TRUE)

  cd <- competing_destinations(pops, D, k = 0.5, a = 1, b = 1, c = 2, delta = 0.5)

  expect_equal(diag(cd), c(0, 0, 0, 0))
  expect_equal(cd, expected, tolerance = 1e-1)
})

test_that("row_normalizer caps row sums at max_rowsum and leaves small rows alone", {
  # Given a matrix with one row summing above the cap (row 1: 0.8) and one below
  #       (row 2: 0.3)
  # When row_normalizer caps rows at 0.5
  # Then row 1 is scaled down to sum exactly 0.5 (proportions preserved) and row 2
  #      is unchanged.
  # Failure would mean over-capping small rows or not capping large ones.
  m <- matrix(c(0.0, 0.6, 0.2,
                0.1, 0.0, 0.2,
                0.3, 0.1, 0.0),
              nrow = 3L, byrow = TRUE)

  capped <- row_normalizer(m, max_rowsum = 0.5)

  expect_equal(sum(capped[1L, ]), 0.5)              # row 1 (was 0.8) capped
  expect_equal(capped[2L, ], m[2L, ])               # row 2 (was 0.3) untouched
  expect_equal(capped[1L, 2L] / capped[1L, 3L], 3)  # proportions preserved (0.6:0.2)
})

test_that("the migration models reject non-symmetric or mis-shaped distances", {
  # Given distance inputs that violate the documented contract
  # When a model is called
  # Then it errors rather than producing nonsense.
  # Failure would risk silent mis-indexing or asymmetric, ill-posed networks.
  asym <- matrix(c(0, 2, 1, 0), nrow = 2L)          # not symmetric
  expect_error(radiation(c(10, 20), asym, k = 1, include_home = FALSE), "symmetric")
  expect_error(gravity(c(10, 20, 30), D, k = 1, a = 1, b = 1, c = 1))  # 3 pops vs 4x4 D
})
