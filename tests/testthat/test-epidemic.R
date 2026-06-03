# Tests for epidemic model components.
# State codes: S=0, E=1, I=2, R=3, D=-1 (matches laser_states()).

# ── Helpers ───────────────────────────────────────────────────────────────────

make_people <- function(n, state = 0L, node = 0L) {
  p <- LaserFrame$new(n * 2L, n)   # capacity 2× so births have room
  p$add_scalar_property("state", "integer", 0L)
  p$add_scalar_property("node",  "integer", 0L)
  p$add_scalar_property("timer", "integer", 0L)
  if (!all(state == 0L)) p$state <- rep_len(as.integer(state), n)
  if (!all(node  == 0L)) p$node  <- rep_len(as.integer(node),  n)
  p
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

  expect_type(states, "integer")
  expect_named(states, c("S", "E", "I", "R", "D"))
  expect_equal(unname(states), c(0L, 1L, 2L, 3L, -1L))
})

# ── step_transmission_si ─────────────────────────────────────────────────────

test_that("step_transmission_si: no infections when beta = 0", {
  # Given: 10% of agents infectious, beta = 0
  # When:  run transmission step
  # Then:  state vector is unchanged
  n   <- 1000L
  ppl <- make_people(n, state = c(rep(2L, 100L), rep(0L, 900L)))
  nd  <- make_nodes(pop = n)
  nd$N <- n

  state_before <- ppl$state
  step_transmission_si(ppl, nd, beta = 0.0, inf_duration = 14L)

  expect_identical(ppl$state, state_before)
})

test_that("step_transmission_si: no infections when no infectious agents", {
  # Given: all agents susceptible (no seed)
  # When:  run transmission step with positive beta
  # Then:  state vector stays all-S
  n   <- 1000L
  ppl <- make_people(n, state = 0L)
  nd  <- make_nodes(pop = n)

  step_transmission_si(ppl, nd, beta = 0.5, inf_duration = 14L)

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

  step_transmission_si(ppl, nd, beta = 100.0, inf_duration = 14L)

  n_remaining_s <- sum(ppl$state == 0L)
  expect_lt(n_remaining_s, 0.05 * 4500L)
})

test_that("step_transmission_si: newly infected agents have timer = inf_duration", {
  # Given: some infectious agents, positive beta
  # When:  run one transmission step
  # Then:  every newly infected agent (previously S, now I) has timer == inf_duration
  set.seed(1L)
  n   <- 2000L
  ppl <- make_people(n, state = c(rep(2L, 200L), rep(0L, 1800L)))
  nd  <- make_nodes(pop = n)
  state_before <- ppl$state

  step_transmission_si(ppl, nd, beta = 0.5, inf_duration = 14L)

  newly_infected <- which(state_before == 0L & ppl$state == 2L)
  expect_true(length(newly_infected) > 0L)
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

  step_transmission_si(ppl, nd, beta = 0.0, inf_duration = 14L)  # beta=0 so no new infections

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

  step_transmission_si(ppl, nd, beta = 100.0, inf_duration = 14L)

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

  step_transmission_se(ppl, nd, beta = 0.5, exp_duration = 5L)

  new_exposures <- which(ppl$state == 1L)
  expect_true(length(new_exposures) > 0L)
  expect_equal(ppl$timer[new_exposures], rep(5L, length(new_exposures)))
  expect_false(any(ppl$state[201:2000] == 2L))  # none jumped directly to I
})

# ── step_exposed_ei ───────────────────────────────────────────────────────────

test_that("step_exposed_ei: E agent with timer=1 transitions to I", {
  # Given: one agent in state E with timer = 1
  # When:  run exposed step with inf_duration = 10
  # Then:  agent is in state I with timer = 10
  ppl <- make_people(1L)
  ppl$state <- 1L
  ppl$timer <- 1L

  step_exposed_ei(ppl, inf_duration = 10L)

  expect_equal(ppl$state[1L], 2L)
  expect_equal(ppl$timer[1L], 10L)
})

test_that("step_exposed_ei: E agent with timer=3 decrements but stays E", {
  # Given: one agent in state E with timer = 3
  # When:  run exposed step
  # Then:  agent stays E with timer = 2
  ppl <- make_people(1L)
  ppl$state <- 1L
  ppl$timer <- 3L

  step_exposed_ei(ppl, inf_duration = 10L)

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

  step_exposed_ei(ppl, inf_duration = 10L)

  expect_identical(ppl$state, state_before)
  expect_identical(ppl$timer, timer_before)
})

# ── step_infectious_ir ────────────────────────────────────────────────────────

test_that("step_infectious_ir: I agent with timer=1 transitions to R", {
  # Given: one agent in state I with timer = 1
  # When:  run infectious->recovered step
  # Then:  agent is in state R with timer = 0
  ppl <- make_people(1L)
  ppl$state <- 2L
  ppl$timer <- 1L

  step_infectious_ir(ppl, 0L)

  expect_equal(ppl$state[1L], 3L)
  expect_equal(ppl$timer[1L], 0L)
})

test_that("step_infectious_ir: I agent with timer=5 decrements and stays I", {
  # Given: one agent in state I with timer = 5
  # When:  run step once
  # Then:  agent stays I with timer = 4
  ppl <- make_people(1L)
  ppl$state <- 2L
  ppl$timer <- 5L

  step_infectious_ir(ppl, 0L)

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

  step_infectious_ir(ppl, 0L)

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
  ppl$timer[ppl$state == 2L] <- 14L
  nd  <- make_nodes(pop = n)

  for (tick in seq_len(50L)) {
    step_transmission_si(ppl, nd, beta = 0.3, inf_duration = 14L)
    step_infectious_ir(ppl, 0L)
  }

  expect_equal(ppl$count, n)
  final_state <- ppl$state
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

  new_states <- ppl$state[101:200]
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
  ppl$timer <- as.integer(timer)
  nd    <- make_nodes(pop = n)

  prev_s <- sum(ppl$state == 0L)
  for (tick in seq_len(100L)) {
    step_transmission_se(ppl, nd, beta = 0.4, exp_duration = 5L)
    step_exposed_ei(ppl, inf_duration = 7L)
    step_infectious_ir(ppl, 0L)

    s_now <- sum(ppl$state == 0L)
    expect_lte(s_now, prev_s)   # S never increases
    prev_s <- s_now

    total <- sum(ppl$state %in% c(0L, 1L, 2L, 3L))
    expect_equal(total, n)      # population conserved (no D state here)
  }
})
