# Tests for carry_forward_states(): an R convenience wrapper that carries a list of
# 2-D census Columns forward one tick (each column tick -> tick+1, via carry_forward)
# and, optionally, totals a list of columns into a running total at tick+1 (e.g. keep
# N = S + I + R current for the FOI denominator). Written given-when-then.

mkvec <- function(col0_by_node, nticks) {
  # A nticks x n_nodes i32 census with tick 0 = col0_by_node, the rest zero.
  n <- length(col0_by_node)
  buf <- allocate_vector("i32", nticks, n)
  buf$set(c(col0_by_node, rep(0L, (nticks - 1L) * n)))
  buf
}

test_that("carry_forward_states copies each census column forward one tick", {
  # Given three 3-tick x 2-node census buffers (S, I, R) seeded only at tick 0
  # When carry_forward_states carries them from tick 0
  # Then each buffer's tick-1 column equals its tick-0 column (S+I+R left to the
  #      dynamics), and tick 2 stays zero. Failure would mean the carry did not
  #      reach every column in the list — losing the incremental census base.
  S <- mkvec(c(10L, 20L), 3L)
  I <- mkvec(c(1L,  2L),  3L)
  R <- mkvec(c(3L,  4L),  3L)

  carry_forward_states(list(S, I, R), 0L)

  expect_equal(S$values()[2L, ], c(10L, 20L))
  expect_equal(I$values()[2L, ], c(1L,  2L))
  expect_equal(R$values()[2L, ], c(3L,  4L))
  expect_equal(S$values()[3L, ], c(0L, 0L))       # tick 2 untouched
})

test_that("carry_forward_states totals the carried columns into N at tick+1", {
  # Given S, I, R census buffers and an N buffer, all 3-tick x 2-node
  # When carry_forward_states carries S/I/R and totals them into N (default summands)
  # Then N's tick-1 column equals S+I+R per node ((10+1+3, 20+2+4) = (14, 26)), while
  #      N's tick-0 column is left as seeded. Failure would mean the population total
  #      used by calc_foi is wrong.
  S <- mkvec(c(10L, 20L), 3L)
  I <- mkvec(c(1L,  2L),  3L)
  R <- mkvec(c(3L,  4L),  3L)
  N <- mkvec(c(14L, 26L), 3L)                       # N_0 = S+I+R at tick 0

  carry_forward_states(list(S, I, R), 0L, total = N)

  expect_equal(N$values()[2L, ], c(14L, 26L))       # 10+1+3, 20+2+4
})

test_that("carry_forward_states totals a subset via the summands argument", {
  # Given a carry list (S, I, R) but a different set of columns to total (just S + I)
  # When carry_forward_states is given an explicit `summands`
  # Then N at tick+1 sums only the summands (10+1, 20+2 = 11, 22), independent of
  #      which columns were carried. Failure would mean summands is ignored and the
  #      total wrongly tracks the carry list.
  S <- mkvec(c(10L, 20L), 2L)
  I <- mkvec(c(1L,  2L),  2L)
  R <- mkvec(c(3L,  4L),  2L)
  N <- mkvec(c(0L,  0L),  2L)

  carry_forward_states(list(S, I, R), 0L, total = N, summands = list(S, I))

  expect_equal(N$values()[2L, ], c(11L, 22L))       # S + I only
})

test_that("carry_forward_states with no total only carries (no error)", {
  # Given a single census buffer and the default total = NULL
  # When carry_forward_states runs
  # Then it carries the column forward and returns invisibly without needing a total.
  # Failure would mean the optional-total branch is mis-guarded.
  S <- mkvec(c(5L, 6L), 2L)

  expect_invisible(carry_forward_states(list(S), 0L))
  expect_equal(S$values()[2L, ], c(5L, 6L))
})
