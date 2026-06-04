# Tests for epidemic model components.
# State codes: S=0, E=1, I=2, R=3, D=-1 (matches laser_states()).
#
# testthat idioms: `test_that("desc", { ... })` blocks with `expect_*` assertions
# (`expect_equal` compares with tolerance, `expect_identical` is exact including
# type/attributes, `expect_true/false`, `expect_lt/gt/lte`). `L` = integer literal.
# The step_* kernels take a Distribution object (e.g. dist_constant(14)).

# ── Helpers ───────────────────────────────────────────────────────────────────

# Run-helpers: `name <- function(...) { ... }` with `=` defaults in the signature
# (state = 0L, node = 0L are supplied only if the caller omits them).
make_people <- function(n, state = 0L, node = 0L) {
  p <- LaserFrame$new(n * 2L, n)   # $new(capacity, count); capacity 2× so births have room
  p$add_scalar_property("state", "integer", 0L)
  p$add_scalar_property("node",  "integer", 0L)
  p$add_scalar_property("timer", "integer", 0L)
  # `if (cond) expr` has no braces when single-statement. `state == 0L` is
  # vectorized, so all(...) collapses to one logical. rep_len(x, n) recycles x to
  # length n; as.integer() coerces (a bare scalar literal default may be a double).
  if (!all(state == 0L)) p$state <- rep_len(as.integer(state), n)
  if (!all(node  == 0L)) p$node  <- rep_len(as.integer(node),  n)
  p   # last expression is the return value
}

make_nodes <- function(n_nodes = 1L, pop = 1000L) {
  nd <- LaserFrame$new(n_nodes, n_nodes)
  nd$add_scalar_property("N", "integer", as.integer(pop))
  nd$add_scalar_property("I", "integer", 0L)
  nd
}

# ── laser_states() ────────────────────────────────────────────────────────────

test_that("laser_states returns a named integer vector with correct codes", {
  # Given: nothing
  # When:  call laser_states()
  # Then:  get a named integer vector with 5 elements and expected values
  states <- laser_states()

  expect_type(states, "integer")                          # asserts the base type
  expect_named(states, c("S", "E", "I", "R", "D"))        # asserts the names() attribute
  # unname() strips names so the comparison is values-only (named vs unnamed
  # vectors are not identical in R).
  expect_equal(unname(states), c(0L, 1L, 2L, 3L, -1L))
})

# ── step_transmission_si ─────────────────────────────────────────────────────

test_that("step_transmission_si: no infections when beta = 0", {
  # Given: 10% of agents infectious, beta = 0
  # When:  run transmission step
  # Then:  state vector is unchanged
  n   <- 1000L
  # c(rep(2L,100L), rep(0L,900L)) builds one length-1000 vector (100 I, then 900 S).
  ppl <- make_people(n, state = c(rep(2L, 100L), rep(0L, 900L)))
  nd  <- make_nodes(pop = n)
  nd$N <- n

  state_before <- ppl$state   # reads the column into an ordinary R vector (a copy)
  step_transmission_si(ppl, nd, beta = 0.0, inf_dist = dist_constant(14))

  # expect_identical is an exact compare (value + type), used here since the step
  # must not have mutated the column at all.
  expect_identical(ppl$state, state_before)
})

test_that("step_transmission_si: no infections when no infectious agents", {
  # Given: all agents susceptible (no seed)
  # When:  run transmission step with positive beta
  # Then:  state vector stays all-S
  n   <- 1000L
  ppl <- make_people(n, state = 0L)
  nd  <- make_nodes(pop = n)

  step_transmission_si(ppl, nd, beta = 0.5, inf_dist = dist_constant(14))

  expect_true(all(ppl$state == 0L))
})

test_that("step_transmission_si: near-certain infection with very high beta", {
  # Given: 10% infectious agents, extremely high beta
  # When:  run transmission step
  # Then:  almost all S agents are infected (≥ 95% of originally-S agents)
  set.seed(42L)
  n   <- 5000L
  ppl <- make_people(n, state = c(rep(2L, 500L), rep(0L, 4500L)))
  nd  <- make_nodes(pop = n)

  step_transmission_si(ppl, nd, beta = 100.0, inf_dist = dist_constant(14))

  n_remaining_s <- sum(ppl$state == 0L)
  expect_lt(n_remaining_s, 0.05 * 4500L)
})

test_that("step_transmission_si: newly infected agents draw timer from inf_dist", {
  # Given: some infectious agents, positive beta
  # When:  run one transmission step with a dist_constant(14) infectious-period distribution
  # Then:  every newly infected agent (previously S, now I) has timer == 14
  set.seed(1L)
  n   <- 2000L
  ppl <- make_people(n, state = c(rep(2L, 200L), rep(0L, 1800L)))
  nd  <- make_nodes(pop = n)
  state_before <- ppl$state

  step_transmission_si(ppl, nd, beta = 0.5, inf_dist = dist_constant(14))

  # which() returns the integer indices where the logical is TRUE; `&` is the
  # vectorized elementwise AND (`&&` would be scalar-only).
  newly_infected <- which(state_before == 0L & ppl$state == 2L)
  expect_true(length(newly_infected) > 0L)
  # ppl$timer[newly_infected] gathers the timer column at those indices.
  expect_true(all(ppl$timer[newly_infected] == 14L))
})

test_that("step_transmission_si: updates nodes$I with infectious count", {
  # Given: 50 I agents all in a single node
  # When:  run transmission step
  # Then:  nodes$I equals 50 (at least after the tally, before any infections)
  #        Note: step first tallies then infects new agents, so nodes$I reflects
  #        the PRE-step I count.
  n   <- 500L
  ppl <- make_people(n, state = c(rep(2L, 50L), rep(0L, 450L)))
  nd  <- make_nodes(pop = n)

  step_transmission_si(ppl, nd, beta = 0.0, inf_dist = dist_constant(14))  # beta=0 so no new infections

  expect_equal(nd$I, 50L)
})

test_that("step_transmission_si: multi-node FOI is node-local", {
  # Given: 1000 agents in 2 nodes (500 each); node 0 has 100 I, node 1 has 0 I
  # When:  run transmission with beta = 100 (near-certain in node 0, impossible in node 1)
  # Then:  virtually all S agents in node 0 get infected; none in node 1
  set.seed(7L)
  n     <- 1000L
  nodes <- c(rep(0L, 500L), rep(1L, 500L))
  state <- c(rep(2L, 100L), rep(0L, 400L),   # node 0: 100 I, 400 S
             rep(0L, 500L))                   # node 1: 500 S

  ppl <- make_people(n, state = state, node = nodes)
  nd  <- LaserFrame$new(2L, 2L)
  nd$add_scalar_property("N", "integer", 500L)
  nd$add_scalar_property("I", "integer", 0L)

  step_transmission_si(ppl, nd, beta = 100.0, inf_dist = dist_constant(14))

  # Logical indexing: ppl$state[nodes == 1L] keeps elements where the mask is TRUE.
  node1_states <- ppl$state[nodes == 1L]
  expect_true(all(node1_states == 0L))   # node 1: nobody infected
  node0_new_inf <- sum(ppl$state[nodes == 0L] == 2L) - 100L  # new infections in node 0
  expect_gt(node0_new_inf, 350L)         # most of node 0's S agents infected
})

# ── step_transmission_se ─────────────────────────────────────────────────────

test_that("step_transmission_se: exposed agents move to E, not I", {
  # Given: some I agents, positive beta, all others S
  # When:  run exposure step
  # Then:  no agent goes directly to I (2); new non-S agents are in state E (1)
  set.seed(3L)
  n   <- 2000L
  ppl <- make_people(n, state = c(rep(2L, 200L), rep(0L, 1800L)))
  nd  <- make_nodes(pop = n)

  step_transmission_se(ppl, nd, beta = 0.5, exp_dist = dist_constant(5))

  new_exposures <- which(ppl$state == 1L)
  expect_true(length(new_exposures) > 0L)
  # compares two vectors elementwise; rep(5L, k) builds the expected timer vector.
  expect_equal(ppl$timer[new_exposures], rep(5L, length(new_exposures)))
  expect_false(any(ppl$state[201:2000] == 2L))  # none jumped directly to I (201:2000 = the originally-S agents)
})

# ── step_exposed_ei ───────────────────────────────────────────────────────────

test_that("step_exposed_ei: E agent with timer=1 transitions to I", {
  # Given: one agent in state E with timer = 1
  # When:  run exposed step with a dist_constant(10) infectious-period distribution
  # Then:  agent is in state I with timer = 10
  ppl <- make_people(1L)
  ppl$state <- 1L   # writes the (length-1) state column; RHS recycles to fill it
  ppl$timer <- 1L

  step_exposed_ei(ppl, inf_dist = dist_constant(10))

  expect_equal(ppl$state[1L], 2L)   # [1L] reads the single agent's value
  expect_equal(ppl$timer[1L], 10L)
})

test_that("step_exposed_ei: E agent with timer=3 decrements but stays E", {
  # Given: one agent in state E with timer = 3
  # When:  run exposed step
  # Then:  agent stays E with timer = 2
  ppl <- make_people(1L)
  ppl$state <- 1L
  ppl$timer <- 3L

  step_exposed_ei(ppl, inf_dist = dist_constant(10))

  expect_equal(ppl$state[1L], 1L)
  expect_equal(ppl$timer[1L], 2L)
})

test_that("step_exposed_ei: non-E agents are unaffected", {
  # Given: agents in states S, I, R, D (none in E)
  # When:  run exposed step
  # Then:  all states and timers unchanged
  ppl <- make_people(4L, state = c(0L, 2L, 3L, -1L))
  ppl$timer <- c(0L, 5L, 3L, 0L)
  state_before <- ppl$state
  timer_before <- ppl$timer

  step_exposed_ei(ppl, inf_dist = dist_constant(10))

  expect_identical(ppl$state, state_before)
  expect_identical(ppl$timer, timer_before)
})

# ── step_infectious_ir ────────────────────────────────────────────────────────

test_that("step_infectious_ir: I agent with timer=1 transitions to R with immunity timer drawn from imm_dist", {
  # Given: one agent in state I with timer = 1
  # When:  run infectious->recovered step with a dist_constant(7) immunity distribution
  # Then:  agent is in state R with its waning timer set to the drawn value (7)
  # Failure would mean the I→R transition does not seed the R waning timer from
  # imm_dist, which step_recovered_rs later counts down in SEIRS.
  ppl <- make_people(1L)
  ppl$state <- 2L
  ppl$timer <- 1L

  step_infectious_ir(ppl, imm_dist = dist_constant(7))

  expect_equal(ppl$state[1L], 3L)
  expect_equal(ppl$timer[1L], 7L)
})

test_that("step_infectious_ir: I agent with timer=5 decrements and stays I", {
  # Given: one agent in state I with timer = 5
  # When:  run step once
  # Then:  agent stays I with timer = 4
  ppl <- make_people(1L)
  ppl$state <- 2L
  ppl$timer <- 5L

  step_infectious_ir(ppl, imm_dist = dist_constant(0))

  expect_equal(ppl$state[1L], 2L)
  expect_equal(ppl$timer[1L], 4L)
})

test_that("step_infectious_ir: non-I agents are unaffected", {
  # Given: agents in states S, E, R, D
  # When:  run step
  # Then:  states and timers unchanged
  ppl <- make_people(4L, state = c(0L, 1L, 3L, -1L))
  ppl$timer <- c(0L, 3L, 5L, 0L)
  state_before <- ppl$state
  timer_before <- ppl$timer

  step_infectious_ir(ppl, imm_dist = dist_constant(0))

  expect_identical(ppl$state, state_before)
  expect_identical(ppl$timer, timer_before)
})

test_that("step_infectious_ir: SIR model conserves population", {
  # Given: mixed S/I population, no births or deaths
  # When:  run 50 ticks of transmission + recovery
  # Then:  total count never changes; I eventually reaches 0
  set.seed(99L)
  n   <- 2000L
  ppl <- make_people(n, state = c(rep(2L, 50L), rep(0L, 1950L)))
  # Masked assignment: R rewrites `ppl$timer[mask] <- 14L` as get-modify-set, so
  # it reads the whole timer column, sets the masked elements, then writes the
  # full vector back via the `$<-` setter. `ppl$state == 2L` is the logical mask.
  ppl$timer[ppl$state == 2L] <- 14L
  nd  <- make_nodes(pop = n)

  for (tick in seq_len(50L)) {
    step_transmission_si(ppl, nd, beta = 0.3, inf_dist = dist_constant(14))
    step_infectious_ir(ppl, imm_dist = dist_constant(0))
  }

  expect_equal(ppl$count, n)   # $count is the live agent count (a frame property)
  final_state <- ppl$state
  # %in% is vectorized set membership (like Python `in`, applied elementwise).
  expect_true(all(final_state %in% c(0L, 2L, 3L)))
  # After 50 ticks the epidemic should have grown (fewer S than initial)
  expect_lt(sum(final_state == 0L), 1950L)
})

# ── step_infectious_is ────────────────────────────────────────────────────────

test_that("step_infectious_is: I agent with timer=1 transitions back to S", {
  # Given: one agent in state I with timer = 1
  # When:  run SIS infectious step
  # Then:  agent is in state S with timer = 0
  ppl <- make_people(1L)
  ppl$state <- 2L
  ppl$timer <- 1L

  step_infectious_is(ppl)

  expect_equal(ppl$state[1L], 0L)
  expect_equal(ppl$timer[1L], 0L)
})

# ── step_recovered_rs ─────────────────────────────────────────────────────────

test_that("step_recovered_rs: R agent with timer=1 becomes S", {
  # Given: one agent in state R with timer = 1
  # When:  run waning-immunity step
  # Then:  agent is in state S with timer = 0
  ppl <- make_people(1L)
  ppl$state <- 3L
  ppl$timer <- 1L

  step_recovered_rs(ppl)

  expect_equal(ppl$state[1L], 0L)
  expect_equal(ppl$timer[1L], 0L)
})

test_that("step_recovered_rs: R agent with timer=3 decrements and stays R", {
  # Given: one agent in state R with timer = 3
  # When:  run step once
  # Then:  agent stays R with timer = 2
  ppl <- make_people(1L)
  ppl$state <- 3L
  ppl$timer <- 3L

  step_recovered_rs(ppl)

  expect_equal(ppl$state[1L], 3L)
  expect_equal(ppl$timer[1L], 2L)
})

# ── step_mortality_cdr ────────────────────────────────────────────────────────

test_that("step_mortality_cdr: cdr=1.0 kills all living agents", {
  # Given: 900 agents in mixed living states (S, I, R); none pre-dead
  # When:  run mortality step with cdr = 1.0
  # Then:  every agent is now dead (state D) and has timer = 0
  # Note:  pre-existing dead agents keep their timers unchanged (tested separately)
  n   <- 900L
  ppl <- make_people(n, state = c(rep(0L, 400L), rep(2L, 300L), rep(3L, 200L)))
  ppl$timer <- rep(5L, n)

  step_mortality_cdr(ppl, cdr = 1.0)

  expect_true(all(ppl$state == -1L))
  expect_true(all(ppl$timer == 0L))
})

test_that("step_mortality_cdr: cdr=0.0 kills nobody", {
  # Given: 1000 living agents
  # When:  run mortality step with cdr = 0
  # Then:  state vector unchanged
  n   <- 1000L
  ppl <- make_people(n, state = 0L)
  state_before <- ppl$state

  step_mortality_cdr(ppl, cdr = 0.0)

  expect_identical(ppl$state, state_before)
})

test_that("step_mortality_cdr: dead agents are not killed again", {
  # Given: all agents already dead (state = D)
  # When:  run mortality step with cdr = 1.0
  # Then:  state and timer vectors are unchanged (D agents are skipped)
  n   <- 100L
  ppl <- make_people(n, state = -1L)
  ppl$timer <- rep(99L, n)
  timer_before <- ppl$timer

  step_mortality_cdr(ppl, cdr = 1.0)

  expect_true(all(ppl$state == -1L))
  expect_identical(ppl$timer, timer_before)
})

test_that("step_mortality_cdr: stochastic mortality rate is approximately correct", {
  # Given: 100 000 S agents, cdr = 0.01
  # When:  run one mortality tick
  # Then:  number of deaths is within 3 standard deviations of expected
  set.seed(17L)
  n   <- 100000L
  ppl <- make_people(n, state = 0L)

  step_mortality_cdr(ppl, cdr = 0.01)

  n_dead     <- sum(ppl$state == -1L)
  expected   <- n * 0.01
  sd_binom   <- sqrt(n * 0.01 * 0.99)
  expect_gt(n_dead, expected - 4 * sd_binom)
  expect_lt(n_dead, expected + 4 * sd_binom)
})

# ── step_births_cbr ───────────────────────────────────────────────────────────

test_that("step_births_cbr: cbr=0 produces no births", {
  # Given: 100 agents, cbr = 0
  # When:  run birth step
  # Then:  count is unchanged
  ppl <- make_people(100L)
  count_before <- ppl$count

  step_births_cbr(ppl, cbr = 0.0)

  expect_equal(ppl$count, count_before)
})

test_that("step_births_cbr: cbr=1 approximately doubles population (capped at capacity)", {
  # Given: 500 agents, capacity 1000
  # When:  run birth step with cbr = 1 (every agent gives birth)
  # Then:  count increases to 1000 (capped by capacity)
  ppl <- LaserFrame$new(1000L, 500L)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 0L)

  step_births_cbr(ppl, cbr = 1.0)

  expect_equal(ppl$count, 1000L)
})

test_that("step_births_cbr: new agents start in state S with timer 0", {
  # Given: 100 I agents, capacity 200, cbr = 1
  # When:  run birth step
  # Then:  the 100 new agents (indices 101:200) all have state = S (0) and timer = 0
  ppl <- LaserFrame$new(200L, 100L)
  ppl$add_scalar_property("state", "integer", 2L)  # default I for existing
  ppl$add_scalar_property("node",  "integer", 0L)
  ppl$add_scalar_property("timer", "integer", 7L)  # existing have timer=7

  step_births_cbr(ppl, cbr = 1.0)

  new_states <- ppl$state[101:200]   # `101:200` indexes the 100 newly-born agents
  new_timers <- ppl$timer[101:200]
  expect_true(all(new_states == 0L))
  expect_true(all(new_timers == 0L))
})

test_that("step_births_cbr: new agents inherit parent node", {
  # Given: 100 agents all in node 5, capacity 200, cbr = 1
  # When:  run birth step
  # Then:  all new agents are also in node 5
  ppl <- LaserFrame$new(200L, 100L)
  ppl$add_scalar_property("state", "integer", 0L)
  ppl$add_scalar_property("node",  "integer", 5L)
  ppl$add_scalar_property("timer", "integer", 0L)

  step_births_cbr(ppl, cbr = 1.0)

  new_nodes <- ppl$node[101:200]
  expect_true(all(new_nodes == 5L))
})

test_that("step_births_cbr: dead agents do not give birth", {
  # Given: all agents in state D (dead), cbr = 1
  # When:  run birth step
  # Then:  count is unchanged
  ppl <- make_people(100L, state = -1L)
  count_before <- ppl$count

  step_births_cbr(ppl, cbr = 1.0)

  expect_equal(ppl$count, count_before)
})

# ── SEIR integration test ─────────────────────────────────────────────────────

test_that("SEIR model: S only decreases, S+E+I+R is conserved over 100 ticks", {
  # Given: 5000 S agents, 10 I agents (seeded), single node SEIR model
  # When:  run 100 ticks with exposure, incubation, recovery
  # Then:  S count never increases; total S+E+I+R stays at 5010
  set.seed(42L)
  n     <- 5010L
  state <- c(rep(2L, 10L), rep(0L, 5000L))
  timer <- c(rep(7L, 10L),  rep(0L, 5000L))  # seeded I agents have 7 ticks left
  ppl   <- make_people(n, state = state)
  ppl$timer <- as.integer(timer)   # as.integer() guarantees the column's integer type
  nd    <- make_nodes(pop = n)

  prev_s <- sum(ppl$state == 0L)
  for (tick in seq_len(100L)) {
    step_transmission_se(ppl, nd, beta = 0.4, exp_dist = dist_constant(5))
    step_exposed_ei(ppl, inf_dist = dist_constant(7))
    step_infectious_ir(ppl, imm_dist = dist_constant(0))

    s_now <- sum(ppl$state == 0L)
    expect_lte(s_now, prev_s)   # S never increases
    prev_s <- s_now

    total <- sum(ppl$state %in% c(0L, 1L, 2L, 3L))
    expect_equal(total, n)      # population conserved (no D state here)
  }
})
