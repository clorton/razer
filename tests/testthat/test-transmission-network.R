# Tests for the spatial `network` coupling added to step_transmission_si /
# step_transmission_se. The network W is an n_nodes x n_nodes matrix where
# W[i, j] is the fraction of node i's force of infection exported to node j, so
# the per-node coupled rate is
#
#     lambda[k] = r[k]*(1 - sum_j W[k,j]) + sum_i r[i]*W[i,k],   r[k] = beta*I[k]/N[k]
#
# and the per-node infection probability is p[k] = 1 - exp(-lambda[k]).
#
# Note: these kernels draw from Rust's thread-local RNG, which `set.seed()` does
# NOT control. Tests therefore assert deterministic facts (p == 0 => exactly zero
# infections) or use large populations so empirical infection fractions sit many
# standard errors from the predicted probability.

states <- laser_states()
S <- states[["S"]]; I <- states[["I"]]

# Build a two-node people frame: node 0 has `i1` infectious among `n1`, node 1 is
# entirely susceptible with `n2` agents.
two_node_people <- function(n1, i1, n2) {
  n <- n1 + n2
  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", S)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)
  ppl$node <- c(rep(0L, n1), rep(1L, n2))
  st <- rep(S, n)
  st[seq_len(i1)] <- I
  ppl$state <- st
  ppl
}

two_node_frame <- function(n1, n2) {
  nd <- LaserFrame$new(2L, 2L)
  nd$add_scalar_property("N", "integer", 0L)
  nd$N <- c(n1, n2)
  nd$add_scalar_property("I", "integer", 0L)
  nd
}

n_infected_in_node <- function(ppl, node) sum(ppl$state[ppl$node == node] == I)

test_that("an unconnected susceptible node receives no force of infection", {
  # Given node 1 (the susceptible node) has no incoming edge from the infectious
  #       node 0 (W[0,1] = 0), and no local infection of its own
  # When  one transmission step runs
  # Then  node 1's force of infection is exactly 0, so zero agents are infected
  #       there — this is deterministic (the p > 0 guard skips them entirely).
  # Failure would mean force leaks across nodes with no connecting edge.
  ppl <- two_node_people(n1 = 100000L, i1 = 20000L, n2 = 200000L)
  nd  <- two_node_frame(100000L, 200000L)
  # node 1 -> node 0 edge only (0.3); nothing points into node 1.
  W <- matrix(c(0, 0,
                0.3, 0), nrow = 2L, byrow = TRUE)
  step_transmission_si(ppl, nd, beta = 0.5, inf_dist = dist_constant(7), network = W)

  expect_equal(n_infected_in_node(ppl, 1L), 0L)
})

test_that("a directed edge leaks infection into an otherwise-empty node", {
  # Given a directed edge node 0 -> node 1 (W[0,1] = 0.5) and infection only in
  #       node 0
  # When  one transmission step runs
  # Then  node 1 acquires infections (its imported force is positive).
  # Failure would mean the network import term (sum_i r[i]*W[i,k]) is not applied.
  ppl <- two_node_people(n1 = 100000L, i1 = 20000L, n2 = 200000L)
  nd  <- two_node_frame(100000L, 200000L)
  W <- matrix(c(0, 0.5,
                0,   0), nrow = 2L, byrow = TRUE)
  step_transmission_si(ppl, nd, beta = 0.5, inf_dist = dist_constant(7), network = W)

  expect_gt(n_infected_in_node(ppl, 1L), 0L)
})

test_that("per-node infection probability matches the redistribution formula", {
  # Given a large two-node population, a directed edge W[0,1] = m, and a non-zero
  #       diagonal on node 0 (which should cancel in the redistribution)
  # When  one transmission step runs
  # Then  the empirical infection fraction in each node matches the predicted
  #       p = 1 - exp(-lambda): node 1 receives r0*m, node 0 retains r0*(1 - m)
  #       regardless of its diagonal entry.
  # Failure would mean the FOI magnitude or diagonal-cancellation is wrong.
  n1 <- 100000L; i1 <- 20000L; n2 <- 500000L
  beta <- 0.3; m <- 0.4
  ppl <- two_node_people(n1, i1, n2)
  nd  <- two_node_frame(n1, n2)
  # Non-zero diagonal on node 0 (0.3) must not change results (self-export cancels).
  W <- matrix(c(0.3, m,
                0,   0), nrow = 2L, byrow = TRUE)
  step_transmission_si(ppl, nd, beta = beta, inf_dist = dist_constant(7), network = W)

  r0 <- beta * i1 / n1
  # Absolute tolerances are generous multiples of the binomial standard error of
  # an empirical fraction (sqrt(p(1-p)/n_S)): ~6 sigma for node 0 (n_S = 80000)
  # and far more for node 1 (n_S = 500000). expect_equal's `tolerance` is relative,
  # so we compare absolute differences explicitly.
  # node 1: all S, imported force only.
  emp_p2 <- n_infected_in_node(ppl, 1L) / n2
  expect_lt(abs(emp_p2 - (1 - exp(-(r0 * m)))), 0.004)
  # node 0: retains (1 - m) of its force; diagonal 0.3 cancels.
  emp_p1 <- (n_infected_in_node(ppl, 0L) - i1) / (n1 - i1)
  expect_lt(abs(emp_p1 - (1 - exp(-(r0 * (1 - m))))), 0.004)
})

test_that("an all-zero network leaves the nodes independent", {
  # Given an all-zero network matrix with infection only in node 0
  # When  a transmission step runs
  # Then  the unconnected node 1 stays empty (deterministically) while node 0
  #       still infects its own susceptibles — i.e. a zero matrix reproduces the
  #       uncoupled, per-node force of infection.
  # Failure would mean a zero matrix unexpectedly couples the nodes.
  ppl <- two_node_people(100000L, 20000L, 200000L)
  nd  <- two_node_frame(100000L, 200000L)
  step_transmission_si(ppl, nd, beta = 0.5, inf_dist = dist_constant(7),
                       network = matrix(0, 2L, 2L))

  expect_equal(n_infected_in_node(ppl, 1L), 0L)   # node 1: no incoming force
  expect_gt(n_infected_in_node(ppl, 0L), 20000L)  # node 0: its own epidemic proceeds
})

test_that("step_transmission_si rejects a mis-shaped network matrix", {
  # Given a network whose dimensions do not match the node count
  # When  the kernel is called
  # Then  it raises an error rather than reading out of bounds.
  # Failure would risk silent memory mis-indexing in the Rust kernel.
  ppl <- two_node_people(1000L, 100L, 1000L)
  nd  <- two_node_frame(1000L, 1000L)
  expect_error(
    step_transmission_si(ppl, nd, beta = 0.3, inf_dist = dist_constant(7),
                         network = matrix(0, 3L, 3L))
  )
})
