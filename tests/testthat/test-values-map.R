# Tests for values_map(): expands a scalar / per-node / per-tick / matrix value
# into a (n_ticks x n_nodes) f64 Column (LASER's ValuesMap), used to build calc_foi()'s beta and
# seasonality grids. $values() reads a 2-D Column back as an n_ticks x n_nodes
# matrix (rows = ticks, columns = nodes). Written given-when-then.

test_that("values_map expands a scalar to a constant grid", {
  # Given a single value and a 3-tick x 2-node grid
  # When values_map runs
  # Then every cell holds that value.
  # Failure would mean the scalar broadcast is wrong.
  g <- values_map(0.5, 3L, 2L)

  expect_s3_class(g, "Column")
  expect_equal(dim(g$values()), c(3L, 2L))
  expect_true(all(g$values() == 0.5))
})

test_that("values_map expands a per-node vector (constant over time)", {
  # Given a length-n_nodes vector
  # When values_map runs
  # Then every tick row equals that per-node vector.
  # Failure would mean per-node values bleed across nodes or vary over time.
  g <- values_map(c(1, 2), 3L, 2L)

  m <- g$values()
  expect_equal(m[1L, ], c(1, 2))
  expect_equal(m[2L, ], c(1, 2))
  expect_equal(m[3L, ], c(1, 2))
})

test_that("values_map expands a per-tick vector (constant over space)", {
  # Given a length-n_ticks vector
  # When values_map runs
  # Then each tick row is that tick's value repeated across all nodes.
  # Failure would mean per-tick values are mis-placed.
  g <- values_map(c(10, 20, 30), 3L, 2L)

  m <- g$values()
  expect_equal(m[1L, ], c(10, 10))
  expect_equal(m[2L, ], c(20, 20))
  expect_equal(m[3L, ], c(30, 30))
})

test_that("values_map uses an n_ticks x n_nodes matrix as-is", {
  # Given a 3 x 2 matrix
  # When values_map runs
  # Then the grid reads back identical to the input matrix.
  # Failure would mean the row/column orientation is transposed.
  input <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3L, ncol = 2L)   # column-major

  g <- values_map(input, 3L, 2L)

  expect_equal(g$values(), input)
})

test_that("values_map rejects a value whose shape matches nothing", {
  # Given a vector that is neither n_nodes nor n_ticks long, and a mis-shaped matrix
  # When values_map is called
  # Then it errors with a clear message rather than silently mis-broadcasting.
  # Failure would risk a wrong-sized parameter grid.
  expect_error(values_map(c(1, 2, 3, 4), 3L, 2L), "scalar")
  expect_error(values_map(matrix(0, 2L, 2L), 3L, 2L), "must be")
})
