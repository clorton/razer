library(razer)

# ── Construction ──────────────────────────────────────────────────────────────

test_that("new() with initial_count = -1 sets count equal to capacity", {
  # Given: a capacity of 10 and the sentinel initial_count
  # When: we construct a frame
  f <- LaserFrame$new(10L, -1L)
  # Then: both count and capacity are 10
  expect_equal(f$count, 10L)
  expect_equal(f$capacity, 10L)
})

test_that("new() with explicit initial_count initialises count correctly", {
  # Given: capacity 20, initial count 5
  # When: we construct a frame
  f <- LaserFrame$new(20L, 5L)
  # Then: count is 5 and capacity is 20
  expect_equal(f$count, 5L)
  expect_equal(f$capacity, 20L)
})

test_that("new() panics on zero or negative capacity", {
  # Given / When / Then: non-positive capacity must be rejected
  expect_error(LaserFrame$new(0L, 0L))
  expect_error(LaserFrame$new(-1L, -1L))
})

test_that("new() panics when initial_count exceeds capacity", {
  # Given: capacity 10
  # When: initial_count is 11
  # Then: an error is raised (overflow check)
  expect_error(LaserFrame$new(10L, 11L))
})

# ── Scalar property registration and access ───────────────────────────────────

test_that("add_scalar_property fills an integer property with the default value", {
  # Given: a frame with 5 active entries
  f <- LaserFrame$new(5L, 5L)
  # When: we add an integer property with default 0
  f$add_scalar_property("age", "integer", 0L)
  # Then: get() returns five zeros
  expect_equal(f$get("age"), rep(0L, 5L))
})

test_that("add_scalar_property fills a real property with the default value", {
  f <- LaserFrame$new(3L, 3L)
  f$add_scalar_property("weight", "real", 1.5)
  expect_equal(f$get("weight"), rep(1.5, 3L))
})

test_that("add_scalar_property fills a logical property with the default value", {
  f <- LaserFrame$new(4L, 4L)
  f$add_scalar_property("alive", "logical", TRUE)
  expect_equal(f$get("alive"), rep(TRUE, 4L))
})

test_that("set() and get() round-trip integer values without loss", {
  # Given: a frame with an integer property
  f <- LaserFrame$new(6L, 6L)
  f$add_scalar_property("id", "integer", 0L)
  vals <- c(10L, 20L, 30L, 40L, 50L, 60L)
  # When: we write those values
  f$set("id", vals)
  # Then: get() returns them unchanged
  expect_equal(f$get("id"), vals)
})

test_that("set() and get() round-trip real values without loss", {
  f <- LaserFrame$new(3L, 3L)
  f$add_scalar_property("x", "real", 0.0)
  vals <- c(1.1, 2.2, 3.3)
  f$set("x", vals)
  expect_equal(f$get("x"), vals)
})

test_that("set() and get() round-trip logical values without loss", {
  f <- LaserFrame$new(4L, 4L)
  f$add_scalar_property("flag", "logical", FALSE)
  vals <- c(TRUE, FALSE, TRUE, FALSE)
  f$set("flag", vals)
  expect_equal(f$get("flag"), vals)
})

test_that("get() returns only the active slice, not the full capacity backing array", {
  # Given: capacity 10 but only 3 entries active
  f <- LaserFrame$new(10L, 3L)
  f$add_scalar_property("age", "integer", 99L)
  # When: we fetch the property
  result <- f$get("age")
  # Then: length is 3, not 10
  expect_length(result, 3L)
})

test_that("get() on a nonexistent property raises an error", {
  f <- LaserFrame$new(5L, 5L)
  expect_error(f$get("nosuchprop"))
})

test_that("adding a property with a duplicate name raises an error", {
  # Given: 'age' already exists as integer
  f <- LaserFrame$new(5L, 5L)
  f$add_scalar_property("age", "integer", 0L)
  # When: we try to add 'age' again as real
  # Then: the duplicate-name guard fires
  expect_error(f$add_scalar_property("age", "real", 0.0))
})

# ── add() — activating new entries ───────────────────────────────────────────

test_that("add() activates entries and returns a 1-based inclusive range", {
  # Given: a frame with 5 active out of 10 capacity
  f <- LaserFrame$new(10L, 5L)
  f$add_scalar_property("age", "integer", 0L)
  # When: we activate 3 more
  range <- f$add(3L)
  # Then: the returned range covers indices 6..8 and count becomes 8
  expect_equal(range, c(6L, 8L))
  expect_equal(f$count, 8L)
})

test_that("add() leaves the existing active slice unchanged", {
  f <- LaserFrame$new(10L, 5L)
  f$add_scalar_property("id", "integer", 0L)
  f$set("id", 1:5)
  f$add(3L)
  # Original five values must be intact
  expect_equal(f$get("id")[1:5], 1:5)
})

test_that("add() new entries carry the property default, not garbage", {
  f <- LaserFrame$new(10L, 5L)
  f$add_scalar_property("score", "integer", -1L)
  range <- f$add(3L)
  result <- f$get("score")
  expect_equal(result[range[1]:range[2]], c(-1L, -1L, -1L))
})

test_that("add() panics when it would exceed capacity", {
  f <- LaserFrame$new(5L, 5L)
  expect_error(f$add(1L))
})

# ── squash() — compaction ─────────────────────────────────────────────────────

test_that("squash() keeps only entries where mask is TRUE", {
  # Given: 6 entries with ids 1..6
  f <- LaserFrame$new(10L, 6L)
  f$add_scalar_property("id", "integer", 0L)
  f$set("id", 1:6)
  # When: we keep only odd-indexed entries
  f$squash(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
  # Then: count is 3 and remaining ids are 1, 3, 5
  expect_equal(f$count, 3L)
  expect_equal(f$get("id"), c(1L, 3L, 5L))
})

test_that("squash() with all-TRUE mask is a no-op", {
  f <- LaserFrame$new(5L, 5L)
  f$add_scalar_property("val", "integer", 0L)
  f$set("val", 1:5)
  f$squash(rep(TRUE, 5L))
  expect_equal(f$count, 5L)
  expect_equal(f$get("val"), 1:5)
})

test_that("squash() with all-FALSE mask empties the frame", {
  f <- LaserFrame$new(5L, 5L)
  f$add_scalar_property("val", "integer", 0L)
  f$squash(rep(FALSE, 5L))
  expect_equal(f$count, 0L)
})

test_that("squash() panics when mask length does not match count", {
  f <- LaserFrame$new(5L, 5L)
  # mask is too short — length mismatch guard
  expect_error(f$squash(c(TRUE, FALSE)))
})

test_that("squash() compacts vector properties column-wise consistently", {
  # Given: 4 entries, a 2-column vector property S, with each column set
  f <- LaserFrame$new(4L, 4L)
  f$add_vector_property("S", 2L, "integer", 0L)
  f$set_col("S", 1L, c(10L, 20L, 30L, 40L))
  f$set_col("S", 2L, c(100L, 200L, 300L, 400L))
  # When: keep entries 1 and 3 (drop 2 and 4)
  f$squash(c(TRUE, FALSE, TRUE, FALSE))
  # Then: each column has the corresponding values
  expect_equal(f$get_col("S", 1L), c(10L, 30L))
  expect_equal(f$get_col("S", 2L), c(100L, 300L))
})

# ── sort_by() — reordering by permutation ────────────────────────────────────

test_that("sort_by() reorders all scalar properties by the given permutation", {
  # Given: a frame with two scalar properties, unsorted
  f <- LaserFrame$new(5L, 5L)
  f$add_scalar_property("a", "integer", 0L)
  f$add_scalar_property("b", "real", 0.0)
  f$set("a", c(5L, 3L, 1L, 4L, 2L))
  f$set("b", c(50.0, 30.0, 10.0, 40.0, 20.0))
  # When: sort ascending by 'a'
  f$sort_by(order(f$get("a")))
  # Then: both properties are in ascending order
  expect_equal(f$get("a"), 1:5)
  expect_equal(f$get("b"), seq(10.0, 50.0, by = 10.0))
})

test_that("sort_by() identity permutation leaves properties unchanged", {
  f <- LaserFrame$new(4L, 4L)
  f$add_scalar_property("v", "integer", 0L)
  f$set("v", c(4L, 2L, 7L, 1L))
  f$sort_by(1:4)
  expect_equal(f$get("v"), c(4L, 2L, 7L, 1L))
})

# ── Vector property access ────────────────────────────────────────────────────

test_that("set_col() and get_col() round-trip values correctly", {
  # Given: a 4-entry frame with a 3-column real vector property
  f <- LaserFrame$new(4L, 4L)
  f$add_vector_property("S", 3L, "real", 0.0)
  vals <- c(10.0, 20.0, 30.0, 40.0)
  # When: we write column 2
  f$set_col("S", 2L, vals)
  # Then: get_col returns those values
  expect_equal(f$get_col("S", 2L), vals)
})

test_that("columns not explicitly set retain the default value", {
  f <- LaserFrame$new(3L, 3L)
  f$add_vector_property("I", 5L, "integer", 0L)
  # Column 3 was never written
  expect_equal(f$get_col("I", 3L), c(0L, 0L, 0L))
})

test_that("get_col() uses 1-based column indexing", {
  f <- LaserFrame$new(2L, 2L)
  f$add_vector_property("X", 3L, "integer", 0L)
  f$set_col("X", 1L, c(11L, 12L))
  f$set_col("X", 3L, c(31L, 32L))
  expect_equal(f$get_col("X", 1L), c(11L, 12L))
  expect_equal(f$get_col("X", 3L), c(31L, 32L))
})

test_that("get_col() panics on out-of-range column index", {
  f <- LaserFrame$new(2L, 2L)
  f$add_vector_property("X", 3L, "integer", 0L)
  # Column 0 is below the 1-based minimum
  expect_error(f$get_col("X", 0L))
  # Column 4 exceeds ncols = 3
  expect_error(f$get_col("X", 4L))
})

test_that("get_matrix() returns the correct dimensions and column values", {
  # Given: 3 active entries, 4-column real vector property
  f <- LaserFrame$new(3L, 3L)
  f$add_vector_property("S", 4L, "real", 0.0)
  for (col in 1:4) {
    f$set_col("S", col, c(col * 1.0, col * 2.0, col * 3.0))
  }
  # When: we fetch the full matrix
  m <- f$get_matrix("S")
  # Then: dimensions are (3 entries × 4 cols) and column values are correct
  expect_equal(dim(m), c(3L, 4L))
  expect_equal(m[, 1], c(1.0, 2.0, 3.0))
  expect_equal(m[, 4], c(4.0, 8.0, 12.0))
})

test_that("get_matrix() reflects the active count after squash()", {
  # Given: 5 entries, vector property, then squash to 3
  f <- LaserFrame$new(5L, 5L)
  f$add_vector_property("S", 2L, "real", 99.0)
  f$add_scalar_property("keep", "logical", TRUE)
  f$set("keep", c(TRUE, FALSE, TRUE, FALSE, TRUE))
  f$squash(f$get("keep"))
  # When: we get the matrix
  m <- f$get_matrix("S")
  # Then: 3 rows, not 5
  expect_equal(nrow(m), 3L)
})

# ── Metadata ──────────────────────────────────────────────────────────────────

test_that("scalar_names() returns alphabetically sorted property names", {
  f <- LaserFrame$new(5L, 5L)
  f$add_scalar_property("zebra", "integer", 0L)
  f$add_scalar_property("apple", "real", 0.0)
  expect_equal(f$scalar_names(), c("apple", "zebra"))
})

test_that("vector_names() returns alphabetically sorted property names", {
  f <- LaserFrame$new(5L, 5L)
  f$add_vector_property("zed", 3L, "real", 0.0)
  f$add_vector_property("aaa", 2L, "integer", 0L)
  expect_equal(f$vector_names(), c("aaa", "zed"))
})

test_that("vector_ncols() returns the column count for a vector property", {
  f <- LaserFrame$new(5L, 5L)
  f$add_vector_property("T", 7L, "real", 0.0)
  expect_equal(f$vector_ncols("T"), 7L)
})

test_that("describe() includes capacity, count, and property names", {
  f <- LaserFrame$new(10L, 5L)
  f$add_scalar_property("age", "integer", 0L)
  f$add_vector_property("S", 3L, "real", 0.0)
  d <- f$describe()
  expect_type(d, "character")
  expect_match(d, "capacity=10")
  expect_match(d, "count=5")
  expect_match(d, "age")
  expect_match(d, "S")
})

# ── Invisible return values (no NULL printed for void methods) ────────────────

test_that("void-returning methods do not auto-print NULL", {
  # Given: a fresh frame
  f <- LaserFrame$new(3L, 3L)
  # When / Then: these calls should return NULL invisibly
  expect_invisible(f$add_scalar_property("x", "integer", 0L))
  expect_invisible(f$add_vector_property("v", 2L, "real", 0.0))
  expect_invisible(f$set("x", 1:3))
  expect_invisible(f$set_col("v", 1L, c(1.0, 2.0, 3.0)))
  expect_invisible(f$squash(rep(TRUE, 3L)))
  expect_invisible(f$sort_by(3:1))
})

# ── Direct property access via $ ──────────────────────────────────────────────

test_that("frame$prop returns the active scalar property slice", {
  # Given: a frame with an integer property set to known values
  f <- LaserFrame$new(5L, 5L)
  f$add_scalar_property("age", "integer", 0L)
  f$set("age", c(10L, 20L, 30L, 40L, 50L))
  # When: we access via $
  result <- f$age
  # Then: same as $get("age")
  expect_equal(result, c(10L, 20L, 30L, 40L, 50L))
})

test_that("frame$prop returns only active entries after squash()", {
  # Given: 6 entries, squash to 3
  f <- LaserFrame$new(6L, 6L)
  f$add_scalar_property("id", "integer", 0L)
  f$set("id", 1:6)
  f$squash(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
  # When: access via $
  # Then: only 3 values
  expect_equal(f$id, c(1L, 3L, 5L))
})

test_that("frame$mat_prop returns the full vector property matrix", {
  # Given: a 3-entry frame with a 2-column vector property
  f <- LaserFrame$new(3L, 3L)
  f$add_vector_property("S", 2L, "integer", 0L)
  f$set_col("S", 1L, c(1L, 2L, 3L))
  f$set_col("S", 2L, c(4L, 5L, 6L))
  # When: access via $
  m <- f$S
  # Then: identical to $get_matrix("S")
  expect_equal(dim(m), c(3L, 2L))
  expect_equal(m[, 1], c(1L, 2L, 3L))
  expect_equal(m[, 2], c(4L, 5L, 6L))
})

test_that("frame$name returns NULL for an unknown name", {
  f <- LaserFrame$new(3L, 3L)
  expect_null(f$no_such_thing)
})

test_that("count and capacity are plain values, not functions", {
  f <- LaserFrame$new(5L, 3L)
  expect_equal(f$count,    3L)
  expect_equal(f$capacity, 5L)
  expect_false(is.function(f$count))
  expect_false(is.function(f$capacity))
})

# ── Direct property assignment via $<- ───────────────────────────────────────

test_that("frame$prop <- values sets a scalar property", {
  # Given: a frame with an integer property initialised to zeros
  f <- LaserFrame$new(4L, 4L)
  f$add_scalar_property("age", "integer", 0L)
  # When: we assign via $<-
  f$age <- c(5L, 10L, 15L, 20L)
  # Then: $get() reflects the new values
  expect_equal(f$get("age"), c(5L, 10L, 15L, 20L))
})

test_that("frame$prop <- values round-trips for real and logical properties", {
  f <- LaserFrame$new(3L, 3L)
  f$add_scalar_property("weight", "real", 0.0)
  f$add_scalar_property("flag", "logical", FALSE)
  f$weight <- c(1.1, 2.2, 3.3)
  f$flag   <- c(TRUE, FALSE, TRUE)
  expect_equal(f$weight, c(1.1, 2.2, 3.3))
  expect_equal(f$flag,   c(TRUE, FALSE, TRUE))
})

test_that("$<- on a vector property raises an informative error", {
  f <- LaserFrame$new(3L, 3L)
  f$add_vector_property("S", 2L, "real", 0.0)
  expect_error(f$S <- matrix(0, 3, 2), "set_col")
})

test_that("$<- on an unknown name raises an error", {
  f <- LaserFrame$new(3L, 3L)
  expect_error(f$no_such_thing <- 42L)
})
