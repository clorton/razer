# Tests for sir_step(): the per-tick SIR recovery kernel. The caller carries the
# census forward (carry_forward) first; sir_step then recovers expired infectious
# agents and applies the I->R delta to column tick+1, plus the recoveries flow.
# Written given-when-then. The parallel-accumulation test checks the per-core tally
# against a serial census of the resulting agent states.

states <- laser_states()                 # c(S=0, E=1, I=2, R=3, D=-1)
S <- states[["S"]]; I <- states[["I"]]; R <- states[["R"]]

mk <- function(dt, x) { col <- allocate_scalar(dt, length(x)); col$set(x); col }

# A zeroed census buffer (n_ticks+1 ticks x n_nodes), with column 0 seeded.
census <- function(nticks, n_nodes, col0 = rep(0L, n_nodes)) {
  buf <- allocate_vector("i32", nticks + 1L, n_nodes)
  buf$set(c(col0, rep(0L, nticks * n_nodes)))   # column 0 is the first n_nodes entries
  buf
}

test_that("sir_step applies the I->R recovery delta to the carried-forward census", {
  # Given 5 agents in 2 nodes (node0: two I with timers 1,2; node1: S, I timer 1, S),
  #       an I/R census seeded at tick 0 and carried forward to tick 1 by the caller
  # When sir_step runs for tick 0
  # Then expired infectious agents recover; the I->R delta is applied to column 1
  #      (I down, R up); the recoveries flow records the per-node counts; and the
  #      resulting census equals a direct census of the agent states.
  # Failure would mean the recovery transition, the delta, or the tally is wrong.
  state  <- mk("u8",  c(I, I, S, I, S))
  timer  <- mk("u8",  c(1, 2, 0, 1, 0))
  nodeid <- mk("u16", c(0, 0, 1, 1, 1))
  Ic <- census(3L, 2L, c(2, 1)); Rc <- census(3L, 2L)
  recoveries <- allocate_vector("i32", 3L, 2L)
  carry_forward(Ic, 0L); carry_forward(Rc, 0L)     # caller carries the census forward

  sir_step(state, timer, nodeid, 5L, Ic, Rc, recoveries, 0L)

  expect_equal(state$values(), c(R, I, S, R, S))   # agents 0 and 3 recovered
  expect_equal(Ic$values()[2L, ], c(1, 0))         # node0: 2-1, node1: 1-1
  expect_equal(Rc$values()[2L, ], c(1, 1))         # one recovery in each node
  expect_equal(recoveries$values()[1L, ], c(1, 1)) # flow recorded for tick 0
  sv <- state$values()
  expect_equal(Ic$values()[2L, ], as.numeric(tabulate(nodeid$values()[sv == I] + 1L, 2L)))
  expect_equal(Rc$values()[2L, ], as.numeric(tabulate(nodeid$values()[sv == R] + 1L, 2L)))
})

test_that("sir_step's parallel recovery tally matches a serial census at scale", {
  # Given 1,000,000 infectious agents spread randomly over 200 nodes, each with a
  #       1-tick timer (so all recover this tick), census carried forward
  # When sir_step runs (work split across cores with per-core node buffers)
  # Then the parallel-accumulated R census and recoveries flow equal a serial
  #      tabulate of the agents by node, with no lost or double-counted events.
  # Failure would expose a race or reduction bug in the parallel accumulation.
  set.seed(1L)
  n_agents <- 1000000L; n_nodes <- 200L
  nid <- sample.int(n_nodes, n_agents, replace = TRUE) - 1L     # 0-based
  state  <- mk("u8",  rep(I, n_agents))
  timer  <- mk("u8",  rep(1L, n_agents))
  nodeid <- mk("u16", nid)
  per_node <- as.numeric(tabulate(nid + 1L, n_nodes))
  Ic <- census(1L, n_nodes, per_node); Rc <- census(1L, n_nodes)
  recoveries <- allocate_vector("i32", 1L, n_nodes)
  carry_forward(Ic, 0L); carry_forward(Rc, 0L)

  sir_step(state, timer, nodeid, n_agents, Ic, Rc, recoveries, 0L)

  expect_true(all(state$values() == R))                 # everyone recovered
  expect_equal(Ic$values()[2L, ], rep(0, n_nodes))      # I drained to zero
  expect_equal(Rc$values()[2L, ], per_node)             # R == per-node agent census
  expect_equal(recoveries$values(), per_node)           # flow == per-node census (1 tick -> vector)
  expect_equal(sum(recoveries$values()), n_agents)
})

test_that("sir_step validates the tick range and census shapes", {
  # Given a tick beyond the census buffer and a mismatched census width
  # When sir_step is called
  # Then it errors rather than writing out of bounds.
  # Failure would risk out-of-bounds memory access.
  state <- mk("u8", c(I, I)); timer <- mk("u8", c(1, 1)); nodeid <- mk("u16", c(0, 1))
  Ic <- census(2L, 2L); Rc <- census(2L, 2L); rec <- allocate_vector("i32", 2L, 2L)
  expect_error(sir_step(state, timer, nodeid, 2L, Ic, Rc, rec, 2L), "out of range")
  Rbad <- census(2L, 3L)
  expect_error(sir_step(state, timer, nodeid, 2L, Ic, Rbad, rec, 0L), "n_nodes")
})
