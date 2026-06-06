# Tests for births(): crude-birth-rate births with newborns entering maternal immunity
# (M). Each living agent births with per-node probability 1-exp(-rate); each newborn
# activates a reserved slot as M with a maternal-waning timer (u16), dob = current tick,
# and a Kaplan-Meier date of death. Returns the new active count, updates the M census at
# tick+1 and the births flow at tick. Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; I <- states[["I"]]; R <- states[["R"]]; M <- states[["M"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }
carried <- function(col0) {
  buf <- allocate_vector("i32", 2L, length(col0)); buf$set(c(col0, col0)); buf
}
rate_grid <- function(r) { g <- allocate_vector("f64", 1L, 1L); g$set(r); g }
# Life table where the whole cohort dies during year 50 (deterministic year of death).
km_year50 <- kaplan_meier_estimator(c(rep(0, 50), rep(100, 51)))

test_that("births activates reserved slots as M with timer, dob and dod set", {
  # Given 4 living agents (S,S,I,R) in node 0 within capacity-10 arrays, and a birth rate
  #       so high every living agent gives birth (p ~ 1), M census carried to tick 1
  # When births runs for tick 0 with a constant maternal timer of 270
  # Then 4 newborns fill slots 5..8 as M with timer 270, dob 0, and a date of death in
  #       [50*365, 51*365) (year-50 life table, born at tick 0); the active count grows to
  #       8; the M census gains 4 and the births flow records 4. Failure would mean a
  #       newborn property, the count bookkeeping, or the census/flow is wrong.
  state  <- mk("u8",  c(S, S, I, R, rep(S, 6L)))
  timer  <- mk("u16", rep(0L, 10L))
  nodeid <- mk("u16", rep(0L, 10L))
  dob    <- mk("i32", rep(0L, 10L))
  dod    <- mk("u32", rep(0L, 10L))
  Mc <- carried(0); births_flow <- allocate_vector("i32", 1L, 1L)

  new_count <- births(state, timer, nodeid, dob, dod, 4L, rate_grid(100),
                      Mc, births_flow, dist_constant(270), km_year50, 0L)

  expect_equal(new_count, 8L)
  expect_equal(state$values(),  c(S, S, I, R, M, M, M, M, S, S))   # slots 5..8 -> M
  expect_equal(timer$values()[5:8], rep(270L, 4L))                # maternal timer
  expect_equal(dob$values()[5:8],   rep(0L, 4L))                  # born this tick
  expect_true(all(dod$values()[5:8] >= 50L * 365L & dod$values()[5:8] < 51L * 365L))
  expect_equal(Mc$values()[2L, ], 4)                              # M census at tick+1
  expect_equal(births_flow$values(), 4)                          # flow at tick
})

test_that("births with a zero rate adds nobody", {
  # Given a zero birth rate
  # When births runs
  # Then the count is unchanged, no slot is activated, and the census/flow stay zero.
  # Failure would mean spurious births at zero rate.
  state  <- mk("u8",  c(S, I, rep(S, 3L))); timer <- mk("u16", rep(0L, 5L))
  nodeid <- mk("u16", rep(0L, 5L)); dob <- mk("i32", rep(0L, 5L)); dod <- mk("u32", rep(0L, 5L))
  Mc <- carried(0); births_flow <- allocate_vector("i32", 1L, 1L)

  new_count <- births(state, timer, nodeid, dob, dod, 2L, rate_grid(0),
                      Mc, births_flow, dist_constant(270), km_year50, 0L)

  expect_equal(new_count, 2L)
  expect_true(all(state$values() == c(S, I, S, S, S)))
  expect_equal(Mc$values()[2L, ], 0); expect_equal(births_flow$values(), 0)
})

test_that("births never exceeds capacity (excess silently dropped)", {
  # Given 4 living agents in capacity-5 arrays (one free slot) and a birth rate so high
  #       all 4 would give birth
  # When births runs
  # Then only one newborn is placed (the single free slot); the count saturates at 5 and
  #       the recorded births match what was actually placed. Failure would risk an
  #       out-of-bounds slot write.
  state  <- mk("u8",  c(S, S, S, S, S)); timer <- mk("u16", rep(0L, 5L))
  nodeid <- mk("u16", rep(0L, 5L)); dob <- mk("i32", rep(0L, 5L)); dod <- mk("u32", rep(0L, 5L))
  Mc <- carried(0); births_flow <- allocate_vector("i32", 1L, 1L)

  new_count <- births(state, timer, nodeid, dob, dod, 4L, rate_grid(100),
                      Mc, births_flow, dist_constant(270), km_year50, 0L)

  expect_equal(new_count, 5L)                  # filled the one free slot, dropped the rest
  expect_equal(state$values()[5L], M)
  expect_equal(births_flow$values(), 1)
})
