# SEIR model: S → E → I → R  (latent period before becoming infectious)
#
# Components (one tick, downstream-first order — see modeling note):
#   1. step_infectious_ir(people, imm_dist)                      — I→R on timer expiry
#   2. step_exposed_ei(people, inf_dist)                         — E→I on timer expiry
#   3. step_transmission_se(people, nodes, beta, exp_dist)      — S→E
#
# State codes: S=0, E=1, I=2, R=3
# The incubation period (exp_duration) delays the infectious wave relative to SIR,
# slowing the epidemic and lowering the peak I count.
#
# testthat idioms: `test_that("desc", { ... })` blocks with `expect_*` assertions.
# See test-SI.R for the shared run-helper conventions (default args, set.seed,
# $new(capacity, count), the column-major dimnamed `traj` matrix, the seq_len()
# tick loop, and sum(state == code) vectorized tallies).

run_seir <- function(n, n_seed = 100L, beta = 0.4, exp_duration = 5L,
                     inf_duration = 7L, nticks = 300L, seed = 42L) {
  set.seed(seed)

  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)

  sv <- rep(0L, n); sv[seq_len(n_seed)] <- 2L; ppl$state <- sv  # seed as I
  tv <- rep(0L, n); tv[seq_len(n_seed)] <- inf_duration; ppl$timer <- tv

  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  traj <- matrix(0L, nrow = nticks + 1L, ncol = 4L,
                 dimnames = list(NULL, c("S", "E", "I", "R")))
  traj[1L, ] <- c(n - n_seed, 0L, n_seed, 0L)

  for (tick in seq_len(nticks)) {
    step_infectious_ir(ppl, imm_dist = dist_constant(0))
    step_exposed_ei(ppl, inf_dist = dist_constant(inf_duration))
    step_transmission_se(ppl, nd, beta = beta, exp_dist = dist_constant(exp_duration), network = matrix(0, 1, 1))
    traj[tick + 1L, ] <- c(sum(ppl$state == 0L), sum(ppl$state == 1L),
                            sum(ppl$state == 2L), sum(ppl$state == 3L))
  }

  traj
}

# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("SEIR: S + E + I + R = N at every tick", {
  # Given: 10 000 agents, 100 seeded I (R0 ≈ 2.8), 300 ticks
  # When:  run SEIR model
  # Then:  all four compartments sum to N at every tick
  traj <- run_seir(n = 10000L)
  expect_true(all(rowSums(traj) == 10000L))
})

test_that("SEIR: S is monotonically non-increasing", {
  # Given: 10 000 agents, R0 ≈ 2.8
  # When:  run 300 ticks
  # Then:  S[t+1] <= S[t] (only S→E transitions affect S)
  traj <- run_seir(n = 10000L)
  expect_true(all(diff(traj[, "S"]) <= 0L))
})

test_that("SEIR: R is monotonically non-decreasing", {
  # Given: 10 000 agents, R0 ≈ 2.8
  # When:  run 300 ticks
  # Then:  R[t+1] >= R[t] (no waning — agents only enter R in SEIR)
  traj <- run_seir(n = 10000L)
  expect_true(all(diff(traj[, "R"]) >= 0L))
})

test_that("SEIR: epidemic burns out (E and I both return to 0)", {
  # Given: 10 000 agents, R0 ≈ 2.8
  # When:  run 300 ticks
  # Then:  E == 0 and I == 0 at the final tick
  # Failure would indicate E→I or I→R transitions are broken.
  traj <- run_seir(n = 10000L)
  expect_equal(traj[nrow(traj), "E"], 0L, ignore_attr = TRUE)
  expect_equal(traj[nrow(traj), "I"], 0L, ignore_attr = TRUE)
})

test_that("SEIR: E peak precedes I peak (incubation delay)", {
  # Given: 10 000 agents, exp_duration = 5, inf_duration = 7
  # When:  run 300 ticks
  # Then:  the tick at which E is maximised is earlier than the tick at which I is maximised
  # This is the defining signature of the latent period: exposed agents accumulate
  # before the infectious wave, so E peaks first.
  traj <- run_seir(n = 10000L)
  # which.max() returns the 1-based index of the first maximum (the peak tick).
  e_peak_tick <- which.max(traj[, "E"])
  i_peak_tick <- which.max(traj[, "I"])
  expect_lt(e_peak_tick, i_peak_tick)   # expect_lt(a, b): assert a < b
})

test_that("SEIR: E is 0 when beta = 0", {
  # Given: 10 000 agents, 100 seeded I (no S→E possible when beta = 0), 300 ticks
  # When:  run SEIR model with beta = 0
  # Then:  E remains 0 throughout; seeded I agents recover without exposing anyone
  traj <- run_seir(n = 10000L, beta = 0.0)
  expect_true(all(traj[, "E"] == 0L))
  expect_true(all(traj[, "S"] == 9900L))
})

# Cross-model comparison: SEIR peak delay vs SIR is tested in test-SIR.R,
# where run_sir() is available (that file defines and sources run_sir).
