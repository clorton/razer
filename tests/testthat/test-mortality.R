# Tests for mortality(): natural-mortality step that retires agents whose date of death
# `dod` (an absolute tick) has arrived. For each living agent with dod <= tick it sets
# state to D and decrements the M/S/E/I/R census it occupied (at column tick+1), adding
# to the per-node deaths flow at column tick. Written given-when-then.
#
# Note: D = -1 is stored in the u8 state Column as 255, so $values() reads a deceased
# agent as 255 (not -1).

states <- laser_states()                       # c(S=0, E=1, I=2, R=3, M=4, D=-1)
S <- states[["S"]]; E <- states[["E"]]; I <- states[["I"]]
R <- states[["R"]]; M <- states[["M"]]
D_U8 <- 255L                                   # how D (-1) reads back from a u8 column

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }
# A 2-column i32 census with columns 0 and 1 both = `col0` (column 1 simulates the
# caller's carry_forward for tick 0).
carried <- function(col0) {
  buf <- allocate_vector("i32", 2L, length(col0)); buf$set(c(col0, col0)); buf
}

test_that("mortality retires due agents and decrements the right compartment", {
  # Given (tick 0, two nodes) node 0 = {I(dod0), S(dod0), R(dod5)} and node 1 =
  #       {M(dod0), S(dod3)}, census carried to tick 1
  # When mortality runs for tick 0
  # Then the agents with dod <= 0 (the I and S in node 0, the M in node 1) become D and
  #      are removed from their compartments at tick 1; the R (dod 5) and the node-1 S
  #      (dod 3) survive; deaths are 2 in node 0 and 1 in node 1. Failure would mean the
  #      due-date test, the compartment decrement, or the death tally is wrong.
  state  <- mk("u8",  c(I, S, R, M, S))
  dod    <- mk("u32", c(0, 0, 5, 0, 3))
  nodeid <- mk("u16", c(0, 0, 0, 1, 1))
  Mc <- carried(c(0, 1)); Sc <- carried(c(1, 1)); Ec <- carried(c(0, 0))
  Ic <- carried(c(1, 0)); Rc <- carried(c(1, 0))
  deaths <- allocate_vector("i32", 1L, 2L)

  mortality(state, dod, nodeid, 5L, Mc, Sc, Ec, Ic, Rc, deaths, 0L)

  expect_equal(state$values(), c(D_U8, D_U8, R, D_U8, S))   # I,S,M retired; R,S survive
  expect_equal(Ic$values()[2L, ], c(0, 0))                  # node0 I: 1 -> 0
  expect_equal(Sc$values()[2L, ], c(0, 1))                  # node0 S: 1 -> 0; node1 S stays
  expect_equal(Mc$values()[2L, ], c(0, 0))                  # node1 M: 1 -> 0
  expect_equal(Rc$values()[2L, ], c(1, 0))                  # R untouched
  expect_equal(deaths$values(), c(2, 1))                    # node0: I+S; node1: M
  # conservation: every death is exactly one compartment decrement
  total_dec <- (Mc$values()[1L,]-Mc$values()[2L,]) + (Sc$values()[1L,]-Sc$values()[2L,]) +
               (Ic$values()[1L,]-Ic$values()[2L,]) + (Rc$values()[1L,]-Rc$values()[2L,])
  expect_equal(deaths$values(), total_dec)
})

test_that("mortality leaves agents whose dod is in the future untouched", {
  # Given agents all with dod > tick
  # When mortality runs
  # Then nobody dies, the census is unchanged, and deaths are zero.
  # Failure would mean premature deaths.
  state  <- mk("u8",  c(S, I, R)); dod <- mk("u32", c(5, 9, 7)); nodeid <- mk("u16", c(0, 0, 0))
  Mc <- carried(0); Sc <- carried(1); Ec <- carried(0); Ic <- carried(1); Rc <- carried(1)
  deaths <- allocate_vector("i32", 1L, 1L)

  mortality(state, dod, nodeid, 3L, Mc, Sc, Ec, Ic, Rc, deaths, 0L)  # all dod > 0

  expect_equal(state$values(), c(S, I, R))
  expect_equal(Sc$values()[2L, ], 1); expect_equal(Ic$values()[2L, ], 1)
  expect_equal(deaths$values(), 0)
})

test_that("mortality does not re-kill or re-count already-deceased agents", {
  # Given an agent already in D (255) alongside a due living agent
  # When mortality runs
  # Then only the living agent is retired and counted; the already-dead agent is left
  #      alone (no double-count). Failure would mean dead agents leak into the deaths
  #      tally or corrupt the census.
  state  <- mk("u8",  c(D_U8, I)); dod <- mk("u32", c(0, 0)); nodeid <- mk("u16", c(0, 0))
  Mc <- carried(0); Sc <- carried(0); Ec <- carried(0); Ic <- carried(1); Rc <- carried(0)
  deaths <- allocate_vector("i32", 1L, 1L)

  mortality(state, dod, nodeid, 2L, Mc, Sc, Ec, Ic, Rc, deaths, 0L)

  expect_equal(state$values(), c(D_U8, D_U8))
  expect_equal(Ic$values()[2L, ], 0)        # the one live I retired
  expect_equal(deaths$values(), 1)          # exactly one death recorded
})

test_that("mortality validates the tick range", {
  # Given a tick whose census column tick+1 does not exist
  # When mortality is called
  # Then it errors rather than writing out of bounds.
  state <- mk("u8", c(I)); dod <- mk("u32", c(0)); nodeid <- mk("u16", c(0))
  Mc <- carried(0); Sc <- carried(0); Ec <- carried(0); Ic <- carried(1); Rc <- carried(0)
  deaths <- allocate_vector("i32", 2L, 1L)
  expect_error(
    mortality(state, dod, nodeid, 1L, Mc, Sc, Ec, Ic, Rc, deaths, 1L),
    "out of range")
})
