# SEIRS model: S → E → I → R → S  (waning immunity)
#
# Components (one tick, downstream-first order — see modeling note):
#   1. step_recovered_rs(people)                                 — R→S on timer expiry
#   2. step_infectious_ir(people, imm_dist)                      — I→R, draws immunity timer
#   3. step_exposed_ei(people, inf_dist)                         — E→I on timer expiry
#   4. step_transmission_se(people, nodes, beta, exp_dist)      — S→E
#
# State codes: S=0, E=1, I=2, R=3
#
# The R→S waning step distinguishes SEIRS from SEIR: recovered agents become
# susceptible again after `imm_duration` ticks, allowing the epidemic to persist
# and potentially produce multiple waves.

run_seirs <- function(n, n_seed = 200L, beta = 0.5, exp_duration = 3L,
                      inf_duration = 7L, imm_duration = 60L,
                      nticks = 500L, seed = 42L) {
  set.seed(seed)

  ppl <- LaserFrame$new(n, n)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)

  sv <- rep(0L, n); sv[seq_len(n_seed)] <- 2L; ppl$state <- sv   # seed as I
  tv <- rep(0L, n); tv[seq_len(n_seed)] <- inf_duration; ppl$timer <- tv

  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  traj <- matrix(0L, nrow = nticks + 1L, ncol = 4L,
                 dimnames = list(NULL, c("S", "E", "I", "R")))
  traj[1L, ] <- c(n - n_seed, 0L, n_seed, 0L)

  for (tick in seq_len(nticks)) {
    step_recovered_rs(ppl)   # R→S waning; counts down the timer drawn from imm_dist by step_infectious_ir
    step_infectious_ir(ppl, imm_dist = dist_constant(imm_duration))
    step_exposed_ei(ppl,    inf_dist = dist_constant(inf_duration))
    step_transmission_se(ppl, nd, beta = beta, exp_dist = dist_constant(exp_duration))
    traj[tick + 1L, ] <- c(sum(ppl$state == 0L), sum(ppl$state == 1L),
                            sum(ppl$state == 2L), sum(ppl$state == 3L))
  }

  traj
}

# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("SEIRS: S + E + I + R = N at every tick", {
  # Given: 20 000 agents, 200 seeded I, 500 ticks
  # When:  run SEIRS model
  # Then:  all four compartments always sum to N
  traj <- run_seirs(n = 20000L)
  expect_true(all(rowSums(traj) == 20000L))
})

test_that("SEIRS: S increases at some point (waning immunity)", {
  # Given: 20 000 agents, imm_duration = 60 ticks
  # When:  run 500 ticks
  # Then:  S[t+1] > S[t] for at least one t, proving R→S transitions occurred
  # In SEIR this never happens; in SEIRS it must once the R pool builds and wanes.
  # Failure indicates step_recovered_rs is not transitioning R agents back to S.
  traj <- run_seirs(n = 20000L)
  expect_true(any(diff(traj[, "S"]) > 0L))
})

test_that("SEIRS: R decreases at some point (waning immunity)", {
  # Given: 20 000 agents, imm_duration = 60 ticks
  # When:  run 500 ticks
  # Then:  R[t+1] < R[t] for at least one t, i.e. R→S waning removes agents from R
  # This is the mirror of the S-increase test: the same R→S transitions that grow S
  # must shrink R.
  traj <- run_seirs(n = 20000L)
  expect_true(any(diff(traj[, "R"]) < 0L))
})

test_that("SEIRS: S does not stay at its post-wave minimum (recovery to S occurs)", {
  # Given: 20 000 agents, imm_duration = 60
  # When:  run 500 ticks
  # Then:  max(S after tick 200) > min(S during ticks 50:200)
  # The epidemic depletes S in the first wave; waning then refills S above the
  # post-wave trough, which would not happen in SEIR.
  traj <- run_seirs(n = 20000L)
  post_wave_trough <- min(traj[50:200, "S"])
  s_after_200      <- max(traj[201:501, "S"])
  expect_gt(s_after_200, post_wave_trough)
})

test_that("SEIRS: epidemic is self-sustaining (I does not permanently reach 0)", {
  # Given: 20 000 agents, beta = 0.5, imm_duration = 20 ticks
  # When:  run 500 ticks
  # Then:  I > 0 at tick 500 (waning refuels transmission, preventing burnout)
  # imm_duration = 20 is short enough that waning overlaps the declining first wave
  # (earliest waning: tick inf_duration + imm_duration = 7 + 20 = 27), so newly
  # susceptible agents are re-exposed while I is still positive.
  # With imm_duration = 60 the first wave ends (~tick 50) before any waning occurs
  # (~tick 110), causing I to reach 0 and the epidemic to die — hence the shorter
  # immunity period here.
  traj <- run_seirs(n = 20000L, imm_duration = 20L)
  expect_gt(traj[501L, "I"], 0L)
})

test_that("SEIRS: waning immunity disabled gives SEIR-like burnout", {
  # Given: SEIRS model with imm_duration >> nticks (immunity never wanes in practice)
  # When:  run 500 ticks
  # Then:  I == 0 at the final tick (no waning → epidemic burns out like SEIR)
  # This confirms that the SEIR→SEIRS distinction is controlled by imm_duration alone.
  traj <- run_seirs(n = 20000L, imm_duration = 9999L)
  expect_equal(traj[501L, "I"], 0L, ignore_attr = TRUE)
})
