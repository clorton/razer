# Tests for the generic single-transition kernels step_timer_expire (untimed destination)
# and step_timer_expire_set (timed destination), the composable building blocks for models
# beyond the named menagerie. Each decrements a u16 timer and, on expiry, moves the agent,
# returning per-node transition counts. Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]; R <- states[["R"]]; M <- states[["M"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }

test_that("step_timer_expire moves from->to (untimed) on timer expiry", {
  # Given two M agents (timers 1 and 3) and an untimed S, in one node
  # When step_timer_expire(M -> S) runs
  # Then the timer-1 agent becomes S (timer 0) and is counted; the timer-3 agent only
  #      decrements; S is untouched. Failure would mean the generic transition is wrong.
  state <- mk("u8", c(M, M, S)); timer <- mk("u16", c(1, 3, 0)); nodeid <- mk("u16", c(0, 0, 0))
  cnt <- step_timer_expire(state, timer, nodeid, 3L, 1L, M, S)

  expect_equal(cnt, 1L)
  expect_equal(state$values(), c(S, M, S))
  expect_equal(timer$values(), c(0, 2, 0))
})

test_that("step_timer_expire_set sets the destination's timer on expiry", {
  # Given two E agents (timers 1 and 4)
  # When step_timer_expire_set(E -> I, constant 9) runs
  # Then the expiring E becomes I with a fresh timer 9 (counted); the other only
  #      decrements. Failure would mean the timed-destination draw is wrong.
  state <- mk("u8", c(E, E)); timer <- mk("u16", c(1, 4)); nodeid <- mk("u16", c(0, 0))
  cnt <- step_timer_expire_set(state, timer, nodeid, 2L, 1L, E, I, dist_constant(9))

  expect_equal(cnt, 1L)
  expect_equal(state$values(), c(I, E))
  expect_equal(timer$values(), c(9, 3))
})

test_that("the generic kernels compose into a hand-built SIRS that conserves population", {
  # Given a single-node SIRS assembled purely from the generic kernels + transmission
  #       (steps first: R->S, I->R, then calc_foi immediately before S->I)
  # When it is run for 80 ticks
  # Then the living compartments sum to N at every tick and an epidemic takes off — i.e.
  #      the building blocks compose into a working model. Failure would mean the generics
  #      cannot be composed correctly.
  set_seed(1)
  n <- 50000L; H <- 80L; nn <- 1L
  state <- allocate_scalar("u8", n); s0 <- rep(S, n); s0[1:50] <- I; state$set(s0)
  timer <- allocate_scalar("u16", n); tm <- rep(0L, n); tm[1:50] <- 7L; timer$set(tm)
  nodeid <- allocate_scalar("u16", n)
  Sc <- allocate_vector("i32", H, nn); Ic <- allocate_vector("i32", H, nn); Rc <- allocate_vector("i32", H, nn)
  N <- allocate_vector("i32", H, nn); foi <- allocate_vector("f64", H - 1L, nn)
  z <- rep(0L, H - 1L); Sc$set(c(n - 50L, z)); Ic$set(c(50L, z)); Rc$set(c(0L, z)); N$set(c(n, z))
  beta <- values_map(0.3, H, nn); season <- values_map(1, H, nn); net <- matrix(0, nn, nn)
  inf <- dist_constant(7); imm <- dist_constant(30)
  for (tick in seq_len(H - 1L)) {
    t <- tick - 1L
    carry_forward_states(list(Sc, Ic, Rc), t, total = N)
    move_count(Rc, Sc, step_timer_expire(state, timer, nodeid, n, nn, R, S), t)          # R->S
    move_count(Ic, Rc, step_timer_expire_set(state, timer, nodeid, n, nn, I, R, imm), t) # I->R (+imm)
    calc_foi(Ic, N, beta, season, net, foi, t)                                            # immediately before transmit
    move_count(Sc, Ic, transmission(state, timer, nodeid, n, foi, t, I, inf), t)         # S->I
  }
  living <- rowSums(Sc$values()) + rowSums(Ic$values()) + rowSums(Rc$values())
  expect_true(all(living == n))
  expect_gt(rowSums(Ic$values())[10L], 50)         # epidemic grew past the seed
})
