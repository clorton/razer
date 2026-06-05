# Tests for distances(): the haversine great-circle distance matrix ported from
# laser-core's `distance` (all-pairs case). Written given-when-then; `L` marks an
# integer literal. distances() is deterministic, so these assert exact-ish values
# (within floating-point / model tolerance), not statistical bounds.

test_that("distances returns a symmetric N x N matrix with a zero diagonal", {
  # Given three well-separated cities (London, Paris, Berlin)
  # When the all-pairs distance matrix is computed
  # Then it is 3 x 3, every diagonal entry (a point to itself) is exactly 0, and
  #      the matrix equals its own transpose (distance is symmetric).
  # Failure would mean the column-major fill or the i<->j indexing is wrong.
  lat <- c(51.5074, 48.8566, 52.5200)
  lon <- c(-0.1278, 2.3522, 13.4050)

  d <- distances(lat, lon)

  expect_equal(dim(d), c(3L, 3L))
  expect_equal(diag(d), c(0, 0, 0))
  expect_equal(d, t(d))
})

test_that("distances matches a known great-circle distance", {
  # Given London and Paris (their reference great-circle separation is ~343.6 km)
  # When the distance matrix is computed
  # Then the off-diagonal entry is within 1 km of the published value, confirming
  #      the haversine formula and the 6371 km Earth radius are wired correctly.
  # Failure would indicate a units error (radians/degrees) or a wrong radius.
  d <- distances(c(51.5074, 48.8566), c(-0.1278, 2.3522))

  expect_lt(abs(d[1L, 2L] - 343.6), 1.0)
})

test_that("distances accepts a single point", {
  # Given a single coordinate
  # When the distance matrix is computed
  # Then the result is the 1 x 1 zero matrix (the point's distance to itself).
  # Failure would mean the degenerate N = 1 case mis-sizes the matrix.
  d <- distances(0, 0)

  expect_equal(dim(d), c(1L, 1L))
  expect_equal(d[1L, 1L], 0)
})

test_that("distances accepts integer coordinate vectors", {
  # Given integer-typed latitude/longitude vectors (R's `L` literals)
  # When the distance matrix is computed
  # Then it succeeds (integers are widened to double in the kernel) and the
  #      equator arc from (0,0) to (0,90) is a quarter of the Earth's
  #      circumference, 6371 * pi/2 ~ 10007.5 km.
  # Failure would mean the integer-slice path is missing or mis-converted.
  d <- distances(c(0L, 0L), c(0L, 90L))

  expect_lt(abs(d[1L, 2L] - 6371 * pi / 2), 1e-6)
})

test_that("distances rejects mismatched-length or out-of-range inputs", {
  # Given coordinate vectors that violate the documented contract
  # When distances() is called
  # Then it raises an error rather than reading out of bounds or returning
  #      nonsense — length mismatch, latitude past the poles, longitude past
  #      the antimeridian.
  # Failure would risk silent mis-indexing or meaningless distances.
  expect_error(distances(c(0, 10), c(0)))                 # length mismatch
  expect_error(distances(c(0, 91), c(0, 0)))              # latitude > 90
  expect_error(distances(c(0, 0), c(0, 200)))             # longitude > 180
})
