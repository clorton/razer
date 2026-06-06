# Tests for the seedable RNG (set_seed / unset_seed). After set_seed(s) every kernel's
# randomness is a deterministic function of s and the order of kernel calls, so a run is
# reproducible (and independent of thread count). Written given-when-then.

# Build an all-susceptible population spread deterministically over `nn` nodes and run one
# transmission step at p = 0.5; returns the per-node infection counts. The node layout
# uses no RNG, so the ONLY randomness is razer's seeded RNG.
N  <- 100000L
NN <- 10L
run_transmission <- function() {
  state  <- allocate_scalar("u8",  N); state$set(rep(0L, N))      # all S
  timer  <- allocate_scalar("u16", N)
  nodeid <- allocate_scalar("u16", N); nodeid$set(rep(0:(NN - 1L), length.out = N))
  foi    <- allocate_vector("f64", 1L, NN); foi$set(rep(log(2), NN))   # p = 1 - exp(-log2) = 0.5
  transmission(state, timer, nodeid, N, foi, 0L, laser_states()[["I"]], dist_constant(5))
}

test_that("set_seed makes a kernel reproducible across runs", {
  # Given the same seed
  # When the same transmission step is run twice
  # Then the per-node infection counts are bit-identical (and a sane fraction infect).
  # Failure would mean the seed does not fully determine the draws.
  set_seed(123); a <- run_transmission()
  set_seed(123); b <- run_transmission()

  expect_identical(a, b)
  expect_gt(sum(a), 0L); expect_lt(sum(a), N)         # ~half infect, not all/none
})

test_that("different seeds give different draws", {
  # Given two different seeds
  # When the step is run under each
  # Then the results differ (overwhelmingly likely across 10 nodes).
  set_seed(123); a <- run_transmission()
  set_seed(999); d <- run_transmission()
  expect_false(identical(a, d))
})

test_that("consecutive seeded calls advance the stream, but the sequence is reproducible", {
  # Given a seed, run two transmission steps in a row (re-susceptible between)
  # When the whole sequence is repeated under the same seed
  # Then each step matches its counterpart, yet the two steps differ from each other (the
  #      per-call counter advanced). Failure would mean either non-reproducibility or that
  #      successive kernel calls reuse the same stream.
  seq_run <- function(seed) {
    set_seed(seed)
    state  <- allocate_scalar("u8",  N); state$set(rep(0L, N))
    timer  <- allocate_scalar("u16", N)
    nodeid <- allocate_scalar("u16", N); nodeid$set(rep(0:(NN - 1L), length.out = N))
    foi    <- allocate_vector("f64", 1L, NN); foi$set(rep(log(2), NN))
    I <- laser_states()[["I"]]
    r1 <- transmission(state, timer, nodeid, N, foi, 0L, I, dist_constant(5))
    state$set(rep(0L, N))                                # reset to all-S for an independent draw
    r2 <- transmission(state, timer, nodeid, N, foi, 0L, I, dist_constant(5))
    list(r1 = r1, r2 = r2)
  }
  x <- seq_run(42); y <- seq_run(42)
  expect_identical(x$r1, y$r1)
  expect_identical(x$r2, y$r2)
  expect_false(identical(x$r1, x$r2))
})

test_that("unset_seed reverts to non-reproducible (entropy) draws", {
  # Given the seed cleared
  # When the step is run twice
  # Then the two runs differ (entropy-seeded). Failure would mean unset_seed did not take.
  unset_seed()
  a <- run_transmission(); b <- run_transmission()
  expect_false(identical(a, b))
})
