# SIR model: S → I → R  (permanent immunity after recovery)
#
# Components (one tick, in order):
#   1. step_transmission_si(people, nodes, beta, inf_dist)      — S→I
#   2. step_infectious_ir(people, imm_dist)                      — I→R on timer expiry
#
# State codes: S=0, I=2, R=3
# R0 = beta * inf_duration (discrete-time approximation)

run_sir <- function(n, n_seed = 100L, beta = 0.3, inf_duration = 14L,
                    nticks = 300L, seed = 42L) {
  set.seed(seed)

  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)

  sv <- rep(0L, n); sv[seq_len(n_seed)] <- 2L; ppl$state <- sv
  tv <- rep(0L, n); tv[seq_len(n_seed)] <- inf_duration; ppl$timer <- tv

  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  traj <- matrix(0L, nrow = nticks + 1L, ncol = 3L,
                 dimnames = list(NULL, c("S", "I", "R")))
  traj[1L, ] <- c(n - n_seed, n_seed, 0L)

  for (tick in seq_len(nticks)) {
    step_transmission_si(ppl, nd, beta = beta, inf_dist = dist_constant(inf_duration))
    step_infectious_ir(ppl, imm_dist = dist_constant(0))
    traj[tick + 1L, ] <- c(sum(ppl$state == 0L),
                            sum(ppl$state == 2L),
                            sum(ppl$state == 3L))
  }

  traj
}

# run_seir is reproduced here (from test-SEIR.R) solely for the cross-model
# comparison test at the bottom of this file; test files do not share definitions.
run_seir <- function(n, n_seed = 100L, beta = 0.4, exp_duration = 5L,
                     inf_duration = 7L, nticks = 300L, seed = 42L) {
  set.seed(seed)
  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)
  sv <- rep(0L, n); sv[seq_len(n_seed)] <- 2L; ppl$state <- sv
  tv <- rep(0L, n); tv[seq_len(n_seed)] <- inf_duration; ppl$timer <- tv
  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)
  traj <- matrix(0L, nrow = nticks + 1L, ncol = 4L,
                 dimnames = list(NULL, c("S", "E", "I", "R")))
  traj[1L, ] <- c(n - n_seed, 0L, n_seed, 0L)
  for (tick in seq_len(nticks)) {
    step_transmission_se(ppl, nd, beta = beta, exp_dist = dist_constant(exp_duration))
    step_exposed_ei(ppl, inf_dist = dist_constant(inf_duration))
    step_infectious_ir(ppl, imm_dist = dist_constant(0))
    traj[tick + 1L, ] <- c(sum(ppl$state == 0L), sum(ppl$state == 1L),
                            sum(ppl$state == 2L), sum(ppl$state == 3L))
  }
  traj
}

# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("SIR: S + I + R = N at every tick", {
  # Given: 10 000 agents, 100 seeded I (R0 ≈ 4.2), 300 ticks
  # When:  run SIR model
  # Then:  S + I + R == N at every tick
  traj <- run_sir(n = 10000L)
  expect_true(all(rowSums(traj) == 10000L))
})

test_that("SIR: S is monotonically non-increasing", {
  # Given: 10 000 agents, R0 ≈ 4.2
  # When:  run 300 ticks
  # Then:  S[t+1] <= S[t] (no route from R back to S in SIR)
  traj <- run_sir(n = 10000L)
  expect_true(all(diff(traj[, "S"]) <= 0L))
})

test_that("SIR: R is monotonically non-decreasing", {
  # Given: 10 000 agents, R0 ≈ 4.2
  # When:  run 300 ticks
  # Then:  R[t+1] >= R[t] (agents only enter R, never leave it in SIR)
  traj <- run_sir(n = 10000L)
  expect_true(all(diff(traj[, "R"]) >= 0L))
})

test_that("SIR: epidemic burns out (I returns to 0)", {
  # Given: 10 000 agents, R0 ≈ 4.2
  # When:  run 300 ticks
  # Then:  I == 0 at the final tick (epidemic cannot persist in a closed SIR system)
  # Failure would indicate recovery (I→R) is not occurring correctly.
  traj <- run_sir(n = 10000L)
  expect_equal(traj[nrow(traj), "I"], 0L, ignore_attr = TRUE)
})

test_that("SIR: some S agents escape infection (herd immunity effect)", {
  # Given: 10 000 agents, R0 ≈ 2.1 (beta=0.15, inf=14; herd threshold ≈ 52%)
  # When:  run until burnout
  # Then:  at least 5% of agents remain S at the end
  # With R0≈4.2 (default beta=0.3) the final attack rate is ~97%, leaving only
  # ~300 survivors; beta=0.15 gives ~80% attack rate and ~2000 survivors.
  # Failure would indicate over-infection or incorrect FOI scaling.
  traj <- run_sir(n = 10000L, beta = 0.15)
  expect_gt(traj[nrow(traj), "S"], 500L)
})

test_that("SIR: I rises then falls (single epidemic wave)", {
  # Given: 10 000 agents, R0 ≈ 4.2
  # When:  run 300 ticks
  # Then:  I peaks strictly above the seed count, then returns to 0
  traj <- run_sir(n = 10000L)
  i_peak <- max(traj[, "I"])
  expect_gt(i_peak, 100L)           # epidemic grew above seed
  expect_equal(traj[nrow(traj), "I"], 0L, ignore_attr = TRUE)  # and then burned out
})

test_that("SIR: beta = 0 produces no new infections", {
  # Given: 10 000 agents, 100 seeded I, beta = 0
  # When:  run 300 ticks (seeded agents recover, no new infections)
  # Then:  S stays at 9900 throughout; I drains into R
  traj <- run_sir(n = 10000L, beta = 0.0)
  expect_true(all(traj[, "S"] == 9900L))
  expect_equal(traj[nrow(traj), "I"], 0L,   ignore_attr = TRUE)  # seeded agents eventually recover
  expect_equal(traj[nrow(traj), "R"], 100L, ignore_attr = TRUE)  # all seed infections moved to R
})

test_that("SIR: recovery timer is respected (I duration ≈ inf_duration)", {
  # Given: single I agent, inf_duration = 10, beta = 0 (no new infections)
  # When:  run 10 ticks
  # Then:  the agent recovers on tick 10 (timer decremented each tick)
  set.seed(1L)
  ppl <- LaserFrame$new(1L, 1L)
  ppl$add_scalar_property("state", "integer", 2L)   # I
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 10L)
  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", 1L)
  nd$add_scalar_property("I", "integer", 0L)

  for (tick in seq_len(9L)) {
    step_transmission_si(ppl, nd, beta = 0.0, inf_dist = dist_constant(10))
    step_infectious_ir(ppl, imm_dist = dist_constant(0))
    expect_equal(ppl$state[1L], 2L, info = paste("tick", tick))
  }
  step_transmission_si(ppl, nd, beta = 0.0, inf_dist = dist_constant(10))
  step_infectious_ir(ppl, imm_dist = dist_constant(0))
  expect_equal(ppl$state[1L], 3L)   # recovered on tick 10
})

test_that("SEIR: incubation delays epidemic peak versus SIR with same R0", {
  # Given: SEIR with beta = 0.4, exp = 5, inf = 7  (R0 ≈ 2.8)
  #        SIR  with beta = 0.2, inf = 14            (R0 ≈ 2.8)
  # When:  run both for 300 ticks
  # Then:  SEIR I-peak tick is later than SIR I-peak tick, due to the E delay
  # Placed here (not test-SEIR.R) so that run_sir() is already defined.
  seir <- run_seir(n = 10000L, beta = 0.4, exp_duration = 5L, inf_duration = 7L, seed = 7L)
  sir  <- run_sir( n = 10000L, beta = 0.2,                    inf_duration = 14L, seed = 7L)
  expect_gt(which.max(seir[, "I"]), which.max(sir[, "I"]))
})
