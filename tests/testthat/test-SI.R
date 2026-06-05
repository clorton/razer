# SI model: S → I  (no recovery — infection is permanent)
#
# Components:
#   step_transmission_si(people, nodes, beta, inf_dist)
#
# State codes: S=0, I=2
# The `inf_dist` argument sets the timer when S→I; because no recovery
# step is called, the timer value has no effect on model dynamics.
#
# testthat idioms used throughout: `test_that("desc", { ... })` registers a test
# block (the `{ }` is an ordinary R expression-block passed as an argument);
# `expect_*` are assertions; `L` suffixes (e.g. `100L`) are integer literals
# (bare `100` is a double in R).

# A run-helper: `name <- function(...) { ... }` defines a closure. Arguments with
# `=` in the signature are defaults (e.g. `beta = 0.3`), supplied only if the
# caller omits them; callers pass them by name (`run_si(n = 10000L)`).
run_si <- function(n, n_seed = 100L, beta = 0.3, nticks = 200L, seed = 42L) {
  set.seed(seed)   # seed the global RNG for reproducibility (state is process-global)

  ppl <- LaserFrame$new(n, n)   # $new(capacity, count) constructs the SoA frame
  # $add_scalar_property(name, type, default) adds one per-agent column.
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)

  sv        <- rep(0L, n)        # rep(x, n) builds a length-n vector of x's
  sv[seq_len(n_seed)] <- 2L      # seq_len(k) == 1:k; index assignment is 1-based
  ppl$state <- sv                # `frame$prop <- v` writes the whole column (see laser_frame.R)

  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", n)
  nd$add_scalar_property("I", "integer", 0L)

  # matrix() is column-major; dimnames=list(rownames, colnames) names the columns
  # so they can be indexed by string (traj[, "S"]). NULL leaves row names unset.
  traj <- matrix(0L, nrow = nticks + 1L, ncol = 2L,
                 dimnames = list(NULL, c("S", "I")))
  traj[1L, ] <- c(n - n_seed, n_seed)   # `m[i, ]` selects a whole row; c() packs the values

  # for (x in vec) is R's foreach; seq_len(nticks) is the 1..nticks index range.
  for (tick in seq_len(nticks)) {
    # step kernels take a Distribution object (dist_constant(n)) for the timer draw.
    step_transmission_si(ppl, nd, beta = beta, inf_dist = dist_constant(n), network = matrix(0, 1, 1))
    # `ppl$state == 0L` is a vectorized elementwise compare → logical vector;
    # sum() of a logical counts the TRUEs (vectorized population tally).
    traj[tick + 1L, ] <- c(sum(ppl$state == 0L), sum(ppl$state == 2L))
  }

  traj   # last expression is the function's return value (no `return` needed)
}

# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("SI: S + I = N at every tick", {
  # Given: 10 000 agents, 100 seeded I, beta = 0.3
  # When:  run SI model for 200 ticks
  # Then:  S + I equals N at every tick (no agents created or destroyed)
  traj <- run_si(n = 10000L)
  # rowSums() collapses each row to its total; `all()` reduces the logical vector
  # to a single TRUE/FALSE. expect_true() asserts that scalar is TRUE.
  expect_true(all(rowSums(traj) == 10000L))
})

test_that("SI: S is monotonically non-increasing", {
  # Given: 10 000 agents, beta = 0.3
  # When:  run 200 ticks
  # Then:  S[t+1] <= S[t] for all t (SI has no route from I back to S)
  traj <- run_si(n = 10000L)
  # traj[, "S"] selects the S column by name; diff() returns lagged differences
  # (x[t+1]-x[t]), so `<= 0` means non-increasing.
  expect_true(all(diff(traj[, "S"]) <= 0L))
})

test_that("SI: I is monotonically non-decreasing", {
  # Given: 10 000 agents, beta = 0.3
  # When:  run 200 ticks
  # Then:  I[t+1] >= I[t] for all t (once infectious, always infectious)
  traj <- run_si(n = 10000L)
  expect_true(all(diff(traj[, "I"]) >= 0L))
})

test_that("SI: beta = 0 produces no new infections", {
  # Given: 10 000 agents, 100 seeded I, beta = 0
  # When:  run 200 ticks
  # Then:  S and I remain unchanged at every tick
  # Failure would indicate the FOI computation ignores beta.
  traj <- run_si(n = 10000L, beta = 0.0)
  expect_true(all(traj[, "S"] == 9900L))
  expect_true(all(traj[, "I"] == 100L))
})

test_that("SI: very high beta infects nearly the whole population", {
  # Given: 10 000 agents, 1% seeded I, beta = 5.0
  # When:  run 100 ticks
  # Then:  fewer than 100 S agents remain (< 1% uninfected)
  # Failure would indicate the FOI→probability conversion or per-agent loop is wrong.
  traj <- run_si(n = 10000L, beta = 5.0, nticks = 100L)
  # expect_lt(a, b) asserts a < b. Row 101 is the final tick (row 1 is t=0).
  expect_lt(traj[101L, "S"], 100L)
})

test_that("SI: epidemic grows in early ticks", {
  # Given: 10 000 agents, 100 seeded I, beta = 0.3
  # When:  run 20 ticks
  # Then:  I at tick 20 > I at tick 1 (epidemic is expanding)
  # Failure would indicate transmission is not occurring.
  traj <- run_si(n = 10000L, nticks = 20L)
  expect_gt(traj[21L, "I"], traj[2L, "I"])   # expect_gt(a, b): assert a > b
})

test_that("SI: nodes$I reflects pre-step infectious count each tick", {
  # Given: single node, 200 I agents initially, beta = 0
  # When:  run 1 tick (no new infections with beta=0)
  # Then:  nodes$I == 200 after the tick (tally is based on pre-FOI state)
  set.seed(1L)
  ppl <- LaserFrame$new(1000L, 1000L)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)
  # `;` separates statements on one line; `1:200` is the inclusive range 1..200.
  sv <- rep(0L, 1000L); sv[1:200] <- 2L; ppl$state <- sv
  nd <- LaserFrame$new(1L, 1L)
  nd$add_scalar_property("N", "integer", 1000L)
  nd$add_scalar_property("I", "integer", 0L)
  step_transmission_si(ppl, nd, beta = 0.0, inf_dist = dist_constant(1), network = matrix(0, 1, 1))
  expect_equal(nd$I, 200L)   # nd$I reads the scalar property back from the frame
})
