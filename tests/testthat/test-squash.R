# Tests for squash(): compact a people environment, reclaiming excluded agents' slots and
# keeping all per-agent Columns row-aligned. Written given-when-then.

states <- laser_states()
S <- states[["S"]]; I <- states[["I"]]; R <- states[["R"]]
D_U8 <- 255L                                   # D (-1) as stored in a u8 state column

test_that("squash drops dead agents and compacts every per-agent column consistently", {
  # Given 5 agents (S, D, I, D, R) with distinct timers and node ids
  # When squash runs with the default (alive) mask
  # Then the 3 survivors (S, I, R) move to the front of EVERY column in order, the count
  #      drops to 3, and the columns stay aligned. Failure would mean the columns desync
  #      or the dead are not reclaimed.
  people <- new.env()
  people$count <- 5L; people$capacity <- 5L
  people$state  <- allocate_scalar("u8",  5L); people$state$set(c(S, D_U8, I, D_U8, R))
  people$timer  <- allocate_scalar("u16", 5L); people$timer$set(c(10, 11, 12, 13, 14))
  people$nodeid <- allocate_scalar("u16", 5L); people$nodeid$set(c(0, 1, 2, 3, 4))

  n <- squash(people)

  expect_equal(n, 3L); expect_equal(people$count, 3L)
  expect_equal(people$state$values()[1:3],  c(S, I, R))      # alive, original order
  expect_equal(people$timer$values()[1:3],  c(10L, 12L, 14L))# aligned with the survivors
  expect_equal(people$nodeid$values()[1:3], c(0L, 2L, 4L))   # aligned
})

test_that("squash accepts an explicit keep mask", {
  # Given an explicit mask
  # When squash runs
  # Then only the flagged agents remain, compacted in order.
  people <- new.env(); people$count <- 4L; people$capacity <- 4L
  people$state <- allocate_scalar("u8",  4L); people$state$set(rep(S, 4L))
  people$x     <- allocate_scalar("i32", 4L); people$x$set(c(10, 20, 30, 40))

  n <- squash(people, keep = c(TRUE, FALSE, TRUE, FALSE))

  expect_equal(n, 2L)
  expect_equal(people$x$values()[1:2], c(10L, 30L))
})

test_that("squash validates the keep length and the column shape", {
  # Given a mask of the wrong length, and a 2-D (non-scalar) column
  # When squash / Column$squash is called
  # Then each errors rather than mis-compacting.
  people <- new.env(); people$count <- 3L
  people$state <- allocate_scalar("u8", 3L); people$state$set(rep(S, 3L))
  expect_error(squash(people, keep = c(TRUE, FALSE)), "length")

  twod <- allocate_vector("i32", 2L, 3L)
  expect_error(twod$squash(c(TRUE, FALSE)), "1-D")
})
