# Tests for allocate_scalar() and the Column it returns: a Rust-owned, dtype-tagged
# 1-D property array exposed to R as an opaque external-pointer handle. The data
# lives in Rust; $values() copies a snapshot back into an R vector for inspection.
# Written given-when-then; `L` marks an integer literal.

all_dtypes <- c("i8", "u8", "i16", "u16", "i32", "u32", "f32", "f64")

test_that("allocate_scalar returns a zero-filled Column for every dtype", {
  # Given each of the eight supported element types
  # When a length-4 column is allocated
  # Then it is a Column handle of the requested dtype and length, and every
  #      element reads back as zero.
  # Failure would mean a dtype branch is mis-wired or not zero-initialized,
  #      corrupting the initial population.
  for (dt in all_dtypes) {
    col <- allocate_scalar(dt, 4L)
    expect_s3_class(col, "Column")
    expect_equal(col$dtype(), dt)
    expect_equal(col$length(), 4L)
    expect_true(all(col$values() == 0))
  }
})

test_that("$values() widens to integer for narrow ints and double for u32/f32/f64", {
  # Given columns whose Rust type does (i32) and does not (u32, f32) fit R's
  #       signed 32-bit integer
  # When the snapshot is read back
  # Then narrow integer types surface as R `integer` and u32/f32/f64 as R
  #      `double` (so u32 values above 2^31-1 are representable).
  # Failure would mean an inspection copy truncates or overflows.
  expect_equal(typeof(allocate_scalar("i8",  1L)$values()), "integer")
  expect_equal(typeof(allocate_scalar("i32", 1L)$values()), "integer")
  expect_equal(typeof(allocate_scalar("u32", 1L)$values()), "double")
  expect_equal(typeof(allocate_scalar("f32", 1L)$values()), "double")

  big <- allocate_scalar("u32", 1L)
  big$set(3000000000)                       # > .Machine$integer.max (2147483647)
  expect_equal(big$values(), 3000000000)
})

test_that("aliases resolve to the canonical dtype", {
  # Given the documented dtype aliases
  # When each is allocated
  # Then it reports the canonical dtype name.
  # Failure would mean the alias table is incomplete.
  expect_equal(allocate_scalar("uint8", 1L)$dtype(), "u8")
  expect_equal(allocate_scalar("raw", 1L)$dtype(), "u8")
  expect_equal(allocate_scalar("integer", 1L)$dtype(), "i32")
  expect_equal(allocate_scalar("double", 1L)$dtype(), "f64")
  expect_equal(allocate_scalar("single", 1L)$dtype(), "f32")
})

test_that("$fill broadcasts a value, truncating toward zero for integer types", {
  # Given a freshly allocated u8 column
  # When $fill is called with 7 and then with a fractional 2.9
  # Then every element becomes 7, then 2 (integer truncation toward zero).
  # Failure would mean fill does not write in place or mishandles the cast.
  col <- allocate_scalar("u8", 3L)

  col$fill(7)
  expect_true(all(col$values() == 7L))
  col$fill(2.9)
  expect_true(all(col$values() == 2L))
})

test_that("$set overwrites from an R vector and round-trips through $values", {
  # Given a length-4 u8 column
  # When $set writes a 4-element vector
  # Then $values returns exactly those values; mutation persists on the same
  #      handle (reference semantics — no copy).
  # Failure would mean set does not write the Rust buffer in place.
  col <- allocate_scalar("u8", 4L)

  col$set(c(1, 2, 3, 250))
  expect_equal(col$values(), c(1L, 2L, 3L, 250L))
})

test_that("a large uint8 allocation has the exact requested length", {
  # Given a multi-million-element request (national-population scale)
  # When a u8 column is allocated
  # Then its length is exact.
  # Failure would indicate a size/index type problem in the allocator.
  n <- 30000000L
  expect_equal(allocate_scalar("u8", n)$length(), n)
})

test_that("allocate_scalar and $set reject contract violations", {
  # Given bad inputs
  # When the allocator or setter is called
  # Then it errors rather than silently mis-allocating or mis-writing.
  # Failure would risk an unusable or wrongly-sized property array.
  expect_error(allocate_scalar("float16", 4L), "unknown dtype")
  expect_error(allocate_scalar("u8", -1L), "non-negative")
  col <- allocate_scalar("u8", 3L)
  expect_error(col$set(c(1, 2)), "length")        # wrong length
})
