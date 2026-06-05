# Tests for carry_forward(): copies column `tick` of a 2-D report Column onto
# column `tick + 1`, seeding the next tick's counts so a dynamics kernel can update
# it in place (count[t+1] = count[t] +/- delta). Written given-when-then.

test_that("carry_forward copies a tick's column onto the next tick", {
  # Given a 3-tick x 2-node counter with tick 0 = (5, 7)
  # When carry_forward runs for tick 0
  # Then tick 1 becomes (5, 7), the source (tick 0) is unchanged, and tick 2 stays 0.
  # Failure would mean the source/destination slice offsets are wrong.
  v <- allocate_vector("i32", 3L, 2L)
  v$set(c(5L, 7L,  0L, 0L,  0L, 0L))

  carry_forward(v, 0L)

  m <- v$values()
  expect_equal(m[1L, ], c(5, 7))   # source unchanged
  expect_equal(m[2L, ], c(5, 7))   # copied to tick 1
  expect_equal(m[3L, ], c(0, 0))   # tick 2 untouched
})

test_that("carry_forward works for any element type", {
  # Given f64 and u8 counters
  # When carried forward
  # Then the next tick equals the source for each — the copy is dtype-agnostic.
  # Failure would mean a Storage variant is missing from the dispatch.
  vf <- allocate_vector("f64", 2L, 3L); vf$set(c(1.5, 2.5, 3.5,  0, 0, 0))
  carry_forward(vf, 0L)
  expect_equal(vf$values()[2L, ], c(1.5, 2.5, 3.5))

  vu <- allocate_vector("u8", 2L, 2L); vu$set(c(10L, 20L,  0L, 0L))
  carry_forward(vu, 0L)
  expect_equal(vu$values()[2L, ], c(10, 20))
})

test_that("carry_forward chains across successive ticks", {
  # Given a 4-tick counter seeded at tick 0
  # When carried forward tick by tick
  # Then the value propagates forward one tick at a time.
  # Failure would mean carry_forward can't be composed across a simulation loop.
  v <- allocate_vector("i32", 4L, 1L); v$set(c(9L, 0L, 0L, 0L))

  carry_forward(v, 0L); carry_forward(v, 1L); carry_forward(v, 2L)

  expect_equal(as.numeric(v$values()), c(9, 9, 9, 9))
})

test_that("carry_forward errors when the destination tick is out of range", {
  # Given a 2-tick counter and a source tick whose successor (tick 2) does not exist
  # When carry_forward is called
  # Then it errors rather than writing past the buffer.
  # Failure would risk an out-of-bounds write.
  v <- allocate_vector("i32", 2L, 2L)
  expect_error(carry_forward(v, 1L), "out of range")
})
