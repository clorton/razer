# Tests for measles_step(): the combined timed-transition kernel — M->S (maternal
# waning), E->I (incubation end), I->R (recovery) — over a uint16 timer, applying census
# deltas at column tick+1. Each agent is processed once (branch on its entry state), so a
# just-arrived I (from E->I) is not also recovered this tick. Written given-when-then.

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]; R <- states[["R"]]; M <- states[["M"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }
carried <- function(col0) {
  buf <- allocate_vector("i32", 2L, length(col0)); buf$set(c(col0, col0)); buf
}

test_that("measles_step transitions M->S, E->I, I->R on timer expiry", {
  # Given one node with M(t1), E(t1), I(t1), plus an untimed S and R, census carried to
  #       tick 1, and a constant infectious period of 5 for the E->I draw
  # When measles_step runs for tick 0
  # Then the M waner becomes S (timer 0), the E case becomes I with a fresh timer 5, the
  #       I case recovers to R (timer 0); S and R are untouched. Census: M 1->0, S 1->2,
  #       E 1->0, I stays 1 (one in, one out), R 1->2. Failure would mean a transition,
  #       the new infectious timer, or a census delta is wrong.
  state  <- mk("u8",  c(M, E, I, S, R))
  timer  <- mk("u16", c(1, 1, 1, 0, 0))
  nodeid <- mk("u16", c(0, 0, 0, 0, 0))
  Mc <- carried(1); Sc <- carried(1); Ec <- carried(1); Ic <- carried(1); Rc <- carried(1)

  measles_step(state, timer, nodeid, 5L, Mc, Sc, Ec, Ic, Rc, dist_constant(5), 0L)

  expect_equal(state$values(), c(S, I, R, S, R))   # M->S, E->I, I->R; S,R unchanged
  expect_equal(timer$values(), c(0, 5, 0, 0, 0))   # E->I drew the constant 5
  expect_equal(Mc$values()[2L, ], 0)
  expect_equal(Sc$values()[2L, ], 2)               # 1 + (M->S)
  expect_equal(Ec$values()[2L, ], 0)
  expect_equal(Ic$values()[2L, ], 1)               # 1 + (E->I) - (I->R)
  expect_equal(Rc$values()[2L, ], 2)               # 1 + (I->R)
  # total living conserved (these transitions move agents, never remove them)
  expect_equal(sum(Mc$values()[2L,], Sc$values()[2L,], Ec$values()[2L,],
                   Ic$values()[2L,], Rc$values()[2L,]), 5)
})

test_that("measles_step only decrements timers that have not yet expired", {
  # Given M(t3), E(t2), I(t5) — none expiring this tick
  # When measles_step runs
  # Then no state changes, each timer is decremented by 1, and the census is unchanged.
  # Failure would mean premature transitions or a missed decrement.
  state  <- mk("u8",  c(M, E, I)); timer <- mk("u16", c(3, 2, 5)); nodeid <- mk("u16", c(0, 0, 0))
  Mc <- carried(1); Sc <- carried(0); Ec <- carried(1); Ic <- carried(1); Rc <- carried(0)

  measles_step(state, timer, nodeid, 3L, Mc, Sc, Ec, Ic, Rc, dist_constant(5), 0L)

  expect_equal(state$values(), c(M, E, I))
  expect_equal(timer$values(), c(2, 1, 4))
  expect_equal(Mc$values()[2L, ], 1); expect_equal(Ec$values()[2L, ], 1)
  expect_equal(Ic$values()[2L, ], 1); expect_equal(Sc$values()[2L, ], 0)
})

test_that("measles_step keeps per-node census in sync at scale", {
  # Given 900,000 agents (300k each in M, E, I, all timer 1) spread over 3 nodes, census
  #       carried to tick 1
  # When measles_step runs (work split across cores)
  # Then per node every M->S, E->I, I->R fires: M empties to S, E empties to I, I (after
  #       losing its originals to R but gaining the E arrivals) equals the per-node E
  #       intake, R gains the original I. The census matches a direct agent census, and
  #       total living is conserved. Failure would expose a race in the reduction.
  set.seed(7L)
  per <- 300000L; nn <- 3L
  nid <- sample.int(nn, 3L * per, replace = TRUE) - 1L
  st  <- c(rep(M, per), rep(E, per), rep(I, per))
  state  <- mk("u8",  st); timer <- mk("u16", rep(1L, 3L * per)); nodeid <- mk("u16", nid)
  cnt <- function(code) as.numeric(tabulate(nid[st == code] + 1L, nn))
  Mc <- carried(cnt(M)); Sc <- carried(rep(0, nn)); Ec <- carried(cnt(E))
  Ic <- carried(cnt(I)); Rc <- carried(rep(0, nn))

  measles_step(state, timer, nodeid, 3L * per, Mc, Sc, Ec, Ic, Rc, dist_constant(6), 0L)

  sv <- state$values()
  expect_equal(Sc$values()[2L, ], as.numeric(tabulate(nid[sv == S] + 1L, nn)))
  expect_equal(Ic$values()[2L, ], as.numeric(tabulate(nid[sv == I] + 1L, nn)))
  expect_equal(Rc$values()[2L, ], as.numeric(tabulate(nid[sv == R] + 1L, nn)))
  expect_true(all(Mc$values()[2L, ] == 0) && all(Ec$values()[2L, ] == 0))
  expect_equal(sum(Sc$values()[2L,], Ic$values()[2L,], Rc$values()[2L,]), 3L * per)
})
