# Tests for calc_foi(): the per-node force-of-infection kernel. Given infectious
# counts, a transmission coefficient beta, and a spatial coupling network, it writes
# the network-redistributed FOI rate into a tick's column of a 2-D report. It reads
# the infectious count from the POST-RECOVERY census column `tick + 1`, so it runs
# after sir_step and excludes agents recovering this tick. Deterministic, so exact
# assertions.
#   r[k]   = beta[k] * infected[k]
#   foi[k] = r[k] * (1 - sum_j W[k,j]) + sum_i r[i] * W[i,k]

# Build an infectious census whose column (tick+1) holds `counts`, so a calc_foi
# call for foi `tick` reads exactly `counts`.
infected_at <- function(tick, counts) {
  n   <- length(counts)
  buf <- allocate_vector("i32", tick + 2L, n)              # columns 0 .. tick+1
  vec <- rep(0L, (tick + 2L) * n)
  vec[((tick + 1L) * n + 1L):((tick + 2L) * n)] <- counts  # place at column tick+1
  buf$set(vec)
  buf
}
# foi report with `nticks` ticks x `n` nodes (f64).
foi_buf <- function(nticks, n) allocate_vector("f64", nticks, n)

test_that("calc_foi with a zero network is the uncoupled local rate beta*I", {
  # Given 2 nodes, post-recovery infectious counts (20, 0), beta 0.5, no coupling
  # When calc_foi runs for tick 0 (reading the census column 1)
  # Then each node's FOI is its own local rate beta*I and no force crosses nodes.
  # Failure would mean the local-rate term or the zero-network path is wrong.
  inf <- infected_at(0L, c(20, 0))
  foi <- foi_buf(1L, 2L)

  calc_foi(inf, beta = 0.5, network = matrix(0, 2L, 2L), foi = foi, tick = 0L)

  expect_equal(foi$values(), c(0.5 * 20, 0))   # 1-tick foi -> plain vector
})

test_that("calc_foi redistributes force across a directed edge", {
  # Given infection only in node 0 and an edge exporting fraction m to node 1
  # When calc_foi runs
  # Then node 0 retains (1 - m) of its local rate and node 1 receives m of it.
  # Failure would mean the export/import redistribution is wrong.
  inf <- infected_at(0L, c(20, 0)); m <- 0.4; r0 <- 0.5 * 20
  W <- matrix(c(0, m,
                0, 0), nrow = 2L, byrow = TRUE)
  foi <- foi_buf(1L, 2L)

  calc_foi(inf, beta = 0.5, network = W, foi = foi, tick = 0L)

  expect_equal(foi$values(), c(r0 * (1 - m), r0 * m))
})

test_that("calc_foi accepts a per-node beta (frequency-dependent via beta = global/N)", {
  # Given populations (100, 50) and per-node beta = beta/N
  # When calc_foi runs with no coupling
  # Then node k's FOI is the frequency-dependent rate beta*I[k]/N[k].
  # Failure would mean the per-node beta branch is mis-applied.
  inf <- infected_at(0L, c(20, 0)); gbeta <- 0.5; N <- c(100, 50)
  foi <- foi_buf(1L, 2L)

  calc_foi(inf, beta = gbeta / N, network = matrix(0, 2L, 2L), foi = foi, tick = 0L)

  expect_equal(foi$values(), c(gbeta * 20 / 100, 0))
})

test_that("calc_foi cancels a network's diagonal (self-export has no net effect)", {
  # Given a non-zero diagonal on node 0 plus the directed edge
  # When calc_foi runs
  # Then the result matches the no-diagonal case (self-export returns and cancels).
  # Failure would mean self-loops leak or double-count force.
  inf <- infected_at(0L, c(20, 0)); m <- 0.4; r0 <- 0.5 * 20
  W <- matrix(c(0.3, m,
                0,   0), nrow = 2L, byrow = TRUE)
  foi <- foi_buf(1L, 2L)

  calc_foi(inf, beta = 0.5, network = W, foi = foi, tick = 0L)

  expect_equal(foi$values(), c(r0 * (1 - m), r0 * m))
})

test_that("calc_foi reads census column tick+1 and writes foi column tick", {
  # Given a census whose column 2 (= tick 1's post-recovery column) holds (10, 0)
  #       and a 3-tick FOI buffer
  # When calc_foi targets tick 1
  # Then foi's tick-1 row holds the FOI from that census column; ticks 0 and 2 stay 0.
  # Failure would mean the read (tick+1) or write (tick) offset is wrong.
  inf <- infected_at(1L, c(10, 0))   # counts placed at column 2
  foi <- foi_buf(3L, 2L)

  calc_foi(inf, beta = 1, network = matrix(0, 2L, 2L), foi = foi, tick = 1L)

  m <- foi$values()
  expect_equal(m[1L, ], c(0, 0))
  expect_equal(m[2L, ], c(10, 0))
  expect_equal(m[3L, ], c(0, 0))
})

test_that("calc_foi validates shapes, tick range, and the foi dtype", {
  # Given inputs that violate the contract
  # When calc_foi is called
  # Then it errors rather than mis-indexing or mis-writing.
  # Failure would risk out-of-bounds access or a wrong-typed write.
  inf <- infected_at(0L, c(1, 2))
  foi <- foi_buf(2L, 2L)
  expect_error(calc_foi(inf, 1, matrix(0, 2L, 2L), foi, 2L), "out of range")        # tick >= nticks
  expect_error(calc_foi(infected_at(0L, c(1, 2, 3)), 1, matrix(0, 2L, 2L), foi, 0L), "n_nodes")  # wrong width
  expect_error(calc_foi(inf, c(1, 2, 3), matrix(0, 2L, 2L), foi, 0L), "length 1 or")  # bad beta length
  expect_error(calc_foi(inf, 1, matrix(0, 3L, 3L), foi, 0L), "matrix")               # wrong network shape
  expect_error(calc_foi(inf, 1, matrix(0, 2L, 2L), allocate_vector("i32", 2L, 2L), 0L), "f64")  # foi not f64
})
