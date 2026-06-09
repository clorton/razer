# Test for laser_states(): the named integer vector of state codes shared
# by the Column-based kernels and the R model scripts. Written given-when-then.

test_that("laser_states returns a named integer vector with correct codes", {
  # Given: nothing
  # When:  call laser_states()
  # Then:  get a named integer vector with 6 elements and expected values
  # Failure would mean the state codes the kernels index by are wrong.
  states <- laser_states()

  expect_type(states, "integer")                          # asserts the base type
  expect_named(states, c("S", "E", "I", "R", "M", "D"))   # asserts the names() attribute
  # unname() strips names so the comparison is values-only (named vs unnamed
  # vectors are not identical in R).
  expect_equal(unname(states), c(0L, 1L, 2L, 3L, 4L, -1L))
})
