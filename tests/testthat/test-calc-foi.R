# Tests for calc_foi(): the per-node force-of-infection kernel. The local rate is
#   r[k] = beta[k] * seasonality[k] * infected[k] / population[k]
# redistributed through the network as foi[k] = r[k]*(1 - sum_j W[k,j]) + sum_i r[i]*W[i,k].
# Index conventions: beta/seasonality are read at the interval column `tick`;
# infected/population (census) at the POST-recovery column `tick+1`; foi written at
# `tick`. Deterministic, so exact assertions.

# Build an f64/i32 Column of `n_cols` x length(vals) with `vals` at column `col_idx`.
at_col <- function(dtype, col_idx, vals, n_cols) {
  n   <- length(vals)
  buf <- allocate_vector(dtype, n_cols, n)
  v   <- rep(0, n_cols * n)
  v[(col_idx * n + 1L):((col_idx + 1L) * n)] <- vals
  buf$set(v)
  buf
}
# Standard tick-0 inputs over 2 nodes: census at column 1, modifiers at column 0.
inf_at1  <- function(vals) at_col("i32", 1L, vals, 2L)
pop_at1  <- function(vals) at_col("i32", 1L, vals, 2L)
grid_at0 <- function(vals) at_col("f64", 0L, vals, 1L)
foi2     <- function() allocate_vector("f64", 1L, 2L)

test_that("calc_foi computes the frequency-dependent local rate with no coupling", {
  # Given infectious (20, 0) of population (100, 50), beta 0.5, seasonality 1, zero net
  # When calc_foi runs for tick 0
  # Then foi[k] = beta * seasonality * I[k] / N[k]: node 0 = 0.5*20/100 = 0.1, node 1 = 0.
  # Failure would mean the rate, the 1/N denominator, or the zero-network path is wrong.
  foi <- foi2()
  calc_foi(inf_at1(c(20, 0)), pop_at1(c(100, 50)), grid_at0(c(0.5, 0.5)),
           grid_at0(c(1, 1)), matrix(0, 2L, 2L), foi, 0L)
  expect_equal(foi$values(), c(0.5 * 20 / 100, 0))   # 1-column foi -> plain vector
})

test_that("calc_foi scales the rate by the seasonality factor", {
  # Given the same inputs but seasonality 2 in node 0
  # When calc_foi runs
  # Then node 0's FOI doubles relative to seasonality 1.
  # Failure would mean the seasonal modifier is ignored or mis-applied.
  foi <- foi2()
  calc_foi(inf_at1(c(20, 0)), pop_at1(c(100, 50)), grid_at0(c(0.5, 0.5)),
           grid_at0(c(2, 1)), matrix(0, 2L, 2L), foi, 0L)
  expect_equal(foi$values(), c(2 * 0.5 * 20 / 100, 0))
})

test_that("calc_foi uses a per-node beta", {
  # Given different beta per node
  # When calc_foi runs with no coupling
  # Then each node's FOI uses its own beta.
  # Failure would mean beta is not read per node.
  foi <- foi2()
  calc_foi(inf_at1(c(20, 10)), pop_at1(c(100, 50)), grid_at0(c(0.5, 0.2)),
           grid_at0(c(1, 1)), matrix(0, 2L, 2L), foi, 0L)
  expect_equal(foi$values(), c(0.5 * 20 / 100, 0.2 * 10 / 50))
})

test_that("calc_foi redistributes force across a directed edge", {
  # Given infection only in node 0 and an edge exporting fraction m to node 1
  # When calc_foi runs
  # Then node 0 retains (1 - m) of its local rate and node 1 receives m of it.
  # Failure would mean the export/import redistribution is wrong.
  m <- 0.4; r0 <- 0.5 * 20 / 100
  W <- matrix(c(0, m, 0, 0), nrow = 2L, byrow = TRUE)
  foi <- foi2()
  calc_foi(inf_at1(c(20, 0)), pop_at1(c(100, 50)), grid_at0(c(0.5, 0.5)),
           grid_at0(c(1, 1)), W, foi, 0L)
  expect_equal(foi$values(), c(r0 * (1 - m), r0 * m))
})

test_that("calc_foi cancels a network's diagonal", {
  # Given a non-zero diagonal on node 0 plus the directed edge
  # When calc_foi runs
  # Then the result matches the no-diagonal case (self-export returns and cancels).
  # Failure would mean self-loops leak or double-count force.
  m <- 0.4; r0 <- 0.5 * 20 / 100
  W <- matrix(c(0.3, m, 0, 0), nrow = 2L, byrow = TRUE)
  foi <- foi2()
  calc_foi(inf_at1(c(20, 0)), pop_at1(c(100, 50)), grid_at0(c(0.5, 0.5)),
           grid_at0(c(1, 1)), W, foi, 0L)
  expect_equal(foi$values(), c(r0 * (1 - m), r0 * m))
})

test_that("calc_foi reads modifiers at tick and census at tick+1, writing foi at tick", {
  # Given a multi-tick setup with beta/seasonality at column 1 and census at column 2
  # When calc_foi targets tick 1
  # Then foi's tick-1 row holds the FOI; ticks 0 and 2 stay zero.
  # Failure would mean a read/write tick offset is wrong.
  infected   <- at_col("i32", 2L, c(10, 0), 3L)   # census at column tick+1 = 2
  population <- at_col("i32", 2L, c(100, 50), 3L)
  beta       <- at_col("f64", 1L, c(0.5, 0.5), 3L) # modifiers at column tick = 1
  season     <- at_col("f64", 1L, c(1, 1), 3L)
  foi        <- allocate_vector("f64", 3L, 2L)

  calc_foi(infected, population, beta, season, matrix(0, 2L, 2L), foi, 1L)

  m <- foi$values()
  expect_equal(m[1L, ], c(0, 0))
  expect_equal(m[2L, ], c(0.5 * 10 / 100, 0))
  expect_equal(m[3L, ], c(0, 0))
})

test_that("calc_foi guards against a zero-population node", {
  # Given a node with zero population
  # When calc_foi runs
  # Then its FOI is 0 rather than NaN/Inf from dividing by zero.
  # Failure would propagate NaN into the transmission probability.
  foi <- foi2()
  calc_foi(inf_at1(c(0, 5)), pop_at1(c(0, 50)), grid_at0(c(0.5, 0.5)),
           grid_at0(c(1, 1)), matrix(0, 2L, 2L), foi, 0L)
  expect_equal(foi$values(), c(0, 0.5 * 5 / 50))
})

test_that("calc_foi validates the tick range, shapes, and network", {
  # Given inputs that violate the contract
  # When calc_foi is called
  # Then it errors rather than mis-indexing.
  # Failure would risk out-of-bounds access.
  foi <- allocate_vector("f64", 2L, 2L)
  ok_inf <- at_col("i32", 1L, c(1, 2), 3L); ok_pop <- at_col("i32", 1L, c(1, 2), 3L)
  ok_b <- grid_at0(c(1, 1)); ok_s <- grid_at0(c(1, 1))
  expect_error(calc_foi(ok_inf, ok_pop, ok_b, ok_s, matrix(0, 2L, 2L), foi, 2L), "out of range")  # tick >= nticks
  expect_error(calc_foi(at_col("i32", 1L, c(1, 2, 3), 3L), ok_pop, ok_b, ok_s,
                        matrix(0, 2L, 2L), foi, 0L), "n_nodes")                                    # wrong width
  expect_error(calc_foi(ok_inf, ok_pop, ok_b, ok_s, matrix(0, 3L, 3L), foi, 0L), "matrix")         # wrong network
})
