# Tests for births(): crude-birth-rate births with newborns entering maternal immunity
# (M). Each living agent births with per-node probability 1-exp(-rate); each newborn
# activates a reserved slot as M with a u16 maternal timer, dob = tick, and a Kaplan-Meier
# date of death. It RETURNS `list(count, born)`: the new active count and per-node births
# (the caller adds them to the M census). Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; I <- states[["I"]]; R <- states[["R"]]; M <- states[["M"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }
rate1 <- function(r) { g <- allocate_vector("f64", 1L, 1L); g$set(r); g }
# Life table where the whole cohort dies during year 50 (deterministic year of death).
km_year50 <- kaplan_meier_estimator(c(rep(0, 50), rep(100, 51)))

test_that("births activates reserved slots as M and returns the grown count + births", {
  # Given 4 living agents (S,S,I,R) in capacity-10 arrays and a birth rate so high every
  #       living agent gives birth (p ~ 1), with a constant maternal timer of 270
  # When births runs for tick 0
  # Then 4 newborns fill slots 5..8 as M with timer 270, dob 0, and a date of death in
  #       [50*365, 51*365); the returned count is 8 and born is 4. Failure would mean a
  #       newborn property, the count bookkeeping, or the returned tally is wrong.
  state  <- mk("u8",  c(S, S, I, R, rep(S, 6L)))
  timer  <- mk("u16", rep(0L, 10L))
  nodeid <- mk("u16", rep(0L, 10L))
  dob    <- mk("i32", rep(0L, 10L))
  dod    <- mk("u32", rep(0L, 10L))

  r <- births(state, timer, nodeid, dob, dod, 4L, 1L, rate1(100),
              dist_constant(270), km_year50, 0L)

  expect_equal(r$count, 8L)
  expect_equal(r$born, 4L)
  expect_equal(state$values(), c(S, S, I, R, M, M, M, M, S, S))
  expect_equal(timer$values()[5:8], rep(270L, 4L))
  expect_equal(dob$values()[5:8],   rep(0L, 4L))
  expect_true(all(dod$values()[5:8] >= 50L * 365L & dod$values()[5:8] < 51L * 365L))
})

test_that("births with a zero rate adds nobody", {
  # Given a zero birth rate
  # When births runs
  # Then the count is unchanged, born is 0, and no slot is activated.
  state  <- mk("u8", c(S, I, rep(S, 3L))); timer <- mk("u16", rep(0L, 5L))
  nodeid <- mk("u16", rep(0L, 5L)); dob <- mk("i32", rep(0L, 5L)); dod <- mk("u32", rep(0L, 5L))

  r <- births(state, timer, nodeid, dob, dod, 2L, 1L, rate1(0),
              dist_constant(270), km_year50, 0L)

  expect_equal(r$count, 2L)
  expect_equal(r$born, 0L)
  expect_true(all(state$values() == c(S, I, S, S, S)))
})

test_that("births never exceeds capacity (excess silently dropped)", {
  # Given 4 living agents in capacity-5 arrays (one free slot) and a birth rate so high
  #       all 4 would give birth
  # When births runs
  # Then only one newborn is placed; the count saturates at 5 and born is 1.
  state  <- mk("u8", rep(S, 5L)); timer <- mk("u16", rep(0L, 5L))
  nodeid <- mk("u16", rep(0L, 5L)); dob <- mk("i32", rep(0L, 5L)); dod <- mk("u32", rep(0L, 5L))

  r <- births(state, timer, nodeid, dob, dod, 4L, 1L, rate1(100),
              dist_constant(270), km_year50, 0L)

  expect_equal(r$count, 5L)
  expect_equal(r$born, 1L)
  expect_equal(state$values()[5L], M)
})
