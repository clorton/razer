# Tests for bincount(): a parallel, NumPy-style histogram that counts how many
# elements of an integer Column fall in each bin, writing into a caller-provided
# counts Column. Written given-when-then; `L` marks an integer literal. bincount is
# deterministic, so these assert exact counts.

make <- function(dtype, x) {
  # Helper: allocate a Column of `dtype` and length `length(x)`, set it to `x`.
  col <- allocate_scalar(dtype, length(x))
  col$set(x)
  col
}

test_that("bincount counts each value into the matching bin", {
  # Given values 0,1,1,2,2,2 and a 3-bin counts buffer
  # When bincount runs
  # Then counts holds c(1, 2, 3) — the number of 0s, 1s and 2s.
  # Failure would mean the per-bin tally or the index mapping is wrong.
  values <- make("u16", c(0, 1, 1, 2, 2, 2))
  counts <- allocate_scalar("i32", 3L)

  bincount(values, 3L, counts)

  expect_equal(counts$values(), c(1L, 2L, 3L))
})

test_that("bincount overwrites the counts buffer on reuse (no stale accumulation)", {
  # Given a counts buffer already holding the result of a prior bincount
  # When bincount is run again into the same buffer
  # Then it yields the same counts, not double — the used bins are zeroed before
  #      the per-thread histograms are accumulated in.
  # Failure would mean the buffer is accumulated onto rather than reset.
  values <- make("u16", c(0, 1, 1, 2, 2, 2))
  counts <- allocate_scalar("i32", 3L)

  bincount(values, 3L, counts)
  bincount(values, 3L, counts)

  expect_equal(counts$values(), c(1L, 2L, 3L))
})

test_that("bincount leaves counts entries at or beyond nbins untouched", {
  # Given a 5-element counts buffer prefilled with 9s and nbins = 3
  # When bincount writes 3 bins
  # Then bins 0..2 hold the counts and bins 3..4 keep their prior 9s.
  # Failure would mean bincount zeroes/writes past nbins.
  values <- make("u16", c(0, 1, 1, 2, 2, 2))
  counts <- make("i32", c(9, 9, 9, 9, 9))

  bincount(values, 3L, counts)

  expect_equal(counts$values(), c(1L, 2L, 3L, 9L, 9L))
})

test_that("bincount works for signed and unsigned value types", {
  # Given the same logical values stored as signed i8 and as unsigned u32
  # When each is binned
  # Then both produce identical counts — the generic kernel handles any integer
  #      width and signedness.
  # Failure would mean a value-type branch is missing or mis-cast.
  expected <- c(1L, 0L, 2L, 1L)   # one 0, no 1, two 2s, one 3
  for (dt in c("i8", "u8", "i16", "u16", "i32", "u32")) {
    counts <- allocate_scalar("i32", 4L)
    bincount(make(dt, c(0, 2, 2, 3)), 4L, counts)
    expect_equal(counts$values(), expected, info = dt)
  }
})

test_that("bincount writes into integer and floating-point counts buffers", {
  # Given counts buffers of various numeric types
  # When bincount runs
  # Then each receives the same counts (cast to its type).
  # Failure would mean a counts-type branch is missing.
  values <- make("u8", c(0, 0, 1))
  for (dt in c("u8", "i16", "u32", "f32", "f64")) {
    counts <- allocate_scalar(dt, 2L)
    bincount(values, 2L, counts)
    expect_equal(counts$values(), c(2, 1), info = dt)
  }
})

test_that("bincount on a large input matches a serial reference (parallel correctness)", {
  # Given a large vector of bin indices in 0..nbins-1
  # When bincount tallies them in parallel
  # Then the result matches R's serial tabulate() exactly — the per-thread
  #      histograms reduce without lost or double-counted increments.
  # Failure would expose a race or a reduction bug in the parallel path.
  set.seed(1L)
  nbins <- 50L
  idx <- sample.int(nbins, size = 2e6, replace = TRUE) - 1L   # 0-based indices
  values <- make("i32", idx)
  counts <- allocate_scalar("i32", nbins)

  bincount(values, nbins, counts)

  expect_equal(counts$values(), tabulate(idx + 1L, nbins = nbins))
})

test_that("bincount validates nbins, buffer length, and value type", {
  # Given contract violations
  # When bincount is called
  # Then it errors rather than mis-indexing or producing garbage.
  # Failure would risk out-of-bounds writes or silently wrong counts.
  values <- make("u16", c(0, 1, 2))
  short  <- allocate_scalar("i32", 2L)
  ok     <- allocate_scalar("i32", 3L)
  expect_error(bincount(values, 3L, short), "at least")     # counts too short
  expect_error(bincount(values, -1L, ok), "non-negative")   # negative nbins
  expect_error(bincount(make("f32", c(0, 1, 2)), 3L, ok), "integer")  # float values
})

# ── bincountw: weighted histogram ────────────────────────────────────────────────

test_that("bincountw sums each element's weight into its bin", {
  # Given values 0,0,1,2,2 with float weights 1.5,2.5,4,1,3 and a 3-bin buffer
  # When bincountw runs
  # Then counts holds the per-bin weight sums c(4, 4, 4): bin 0 = 1.5+2.5, bin 1 =
  #      4, bin 2 = 1+3.
  # Failure would mean the weighted accumulation or the index mapping is wrong.
  values  <- make("u16", c(0, 0, 1, 2, 2))
  weights <- make("f64", c(1.5, 2.5, 4, 1, 3))
  counts  <- allocate_scalar("f64", 3L)

  bincountw(values, weights, 3L, counts)

  expect_equal(counts$values(), c(4, 4, 4))
})

test_that("bincountw handles signed, unsigned, and floating-point weights", {
  # Given the same values binned with i32 (incl. a negative), u32, and f32 weights
  # When each is summed
  # Then the per-bin sums reflect the weight type: the signed case can go negative,
  #      and all widen to f64 for accumulation.
  # Failure would mean a weight-type branch is missing or mis-cast.
  values <- make("u8", c(0, 1, 1, 2))

  ci <- allocate_scalar("f64", 3L)
  bincountw(values, make("i32", c(10, -1, 5, 7)), 3L, ci)
  expect_equal(ci$values(), c(10, 4, 7))           # bin 1 = -1 + 5

  cu <- allocate_scalar("f64", 3L)
  bincountw(values, make("u32", c(10, 1, 5, 7)), 3L, cu)
  expect_equal(cu$values(), c(10, 6, 7))

  cf <- allocate_scalar("f64", 3L)
  bincountw(values, make("f32", c(0.5, 1, 0.25, 2)), 3L, cf)
  expect_equal(cf$values(), c(0.5, 1.25, 2))
})

test_that("bincountw truncates toward zero when counts is integer-typed", {
  # Given fractional weights summed into an integer counts buffer
  # When bincountw runs
  # Then each bin's float sum is truncated toward zero on write (3.7 -> 3).
  # Failure would mean the integer cast on write is wrong.
  values  <- make("u8", c(0, 0))
  weights <- make("f64", c(1.9, 1.8))   # sum 3.7
  counts  <- allocate_scalar("i32", 1L)

  bincountw(values, weights, 1L, counts)

  expect_equal(counts$values(), 3L)
})

test_that("bincountw overwrites the counts buffer on reuse and preserves the tail", {
  # Given a reused 4-element counts buffer prefilled with 9s and nbins = 3
  # When bincountw runs twice
  # Then bins 0..2 hold the (same) weighted sums and bin 3 keeps its prior 9.
  # Failure would mean accumulation onto stale data or writing past nbins.
  values  <- make("u16", c(0, 1, 2, 2))
  weights <- make("f64", c(2, 3, 1, 1))
  counts  <- make("f64", c(9, 9, 9, 9))

  bincountw(values, weights, 3L, counts)
  bincountw(values, weights, 3L, counts)

  expect_equal(counts$values(), c(2, 3, 2, 9))
})

test_that("bincountw on a large input matches a serial weighted tally", {
  # Given 2 million (bin index, weight) pairs
  # When bincountw sums the weights per bin in parallel
  # Then the result matches a serial tapply()-based sum to floating tolerance —
  #      the private per-thread accumulators reduce without lost contributions.
  # Failure would expose a race or reduction bug in the weighted path.
  set.seed(2L)
  nbins <- 40L
  idx   <- sample.int(nbins, size = 2e6, replace = TRUE) - 1L
  wt    <- runif(2e6)
  counts <- allocate_scalar("f64", nbins)

  bincountw(make("i32", idx), make("f64", wt), nbins, counts)

  ref <- as.numeric(tapply(wt, factor(idx, levels = 0:(nbins - 1L)), sum))
  ref[is.na(ref)] <- 0
  expect_equal(counts$values(), ref, tolerance = 1e-9)
})

test_that("bincount writes into the requested slot of a 2-D counts buffer", {
  # Given a 3-tick x 2-node report and values over 2 bins
  # When bincount targets slot (tick) 1
  # Then tick 1's row holds the counts and ticks 0 and 2 stay zero; calling with
  #      the default slot writes tick 0 instead.
  # Failure would mean the slot offset into the counts buffer is wrong.
  values <- make("u16", c(0, 1, 1, 1))          # bin 0: one, bin 1: three
  report <- allocate_vector("i32", 3L, 2L)      # 3 ticks x 2 nodes

  bincount(values, 2L, report, 1L)
  m <- report$values()
  expect_equal(m[1L, ], c(0, 0))                # tick 0 untouched
  expect_equal(m[2L, ], c(1, 3))                # tick 1 written
  expect_equal(m[3L, ], c(0, 0))                # tick 2 untouched

  report0 <- allocate_vector("i32", 3L, 2L)
  bincount(values, 2L, report0)                 # slot defaults to 0
  expect_equal(report0$values()[1L, ], c(1, 3))
})

test_that("bincountw writes into the requested slot of a 2-D counts buffer", {
  # Given a 3-tick x 2-node report, values over 2 bins, and float weights
  # When bincountw targets slot (tick) 2
  # Then tick 2's row holds the per-bin weight sums and the other ticks stay zero.
  # Failure would mean the weighted slot offset is wrong.
  values  <- make("u16", c(0, 1, 1, 1))
  weights <- make("f64", c(2, 3, 1, 1))         # bin 0: 2, bin 1: 3+1+1 = 5
  report  <- allocate_vector("f64", 3L, 2L)

  bincountw(values, weights, 2L, report, 2L)
  m <- report$values()
  expect_equal(m[1L, ], c(0, 0))
  expect_equal(m[2L, ], c(0, 0))
  expect_equal(m[3L, ], c(2, 5))
})

test_that("bincount rejects an out-of-range slot", {
  # Given a 2-tick report and a slot index of 2 (valid slots are 0 and 1)
  # When bincount is called
  # Then it errors rather than writing past the buffer.
  # Failure would risk an out-of-bounds write.
  values <- make("u16", c(0, 1))
  report <- allocate_vector("i32", 2L, 2L)
  expect_error(bincount(values, 2L, report, 2L), "out of range")
})

test_that("bincountw validates lengths, nbins, buffer size, and value type", {
  # Given contract violations
  # When bincountw is called
  # Then each raises an error rather than mis-indexing or summing garbage.
  # Failure would risk out-of-bounds writes or silently wrong sums.
  values  <- make("u16", c(0, 1, 2))
  weights <- make("f64", c(1, 1, 1))
  ok      <- allocate_scalar("f64", 3L)
  expect_error(bincountw(values, make("f64", c(1, 1)), 3L, ok), "same length")
  expect_error(bincountw(values, weights, 3L, allocate_scalar("f64", 2L)), "at least")
  expect_error(bincountw(values, weights, -1L, ok), "non-negative")
  expect_error(bincountw(make("f32", c(0, 1, 2)), weights, 3L, ok), "integer")
})

test_that("bincount tallies agent node ids into a per-node census slot (column-model use)", {
  # Given a per-agent u16 `nodeid` Column and a 2-D (ticks x nodes) i32 census buffer
  # When bincount writes the per-node agent counts into one tick's slot
  # Then that slot holds the node histogram (matching tabulate) and the other slots are
  #      untouched. This is the Column-model aggregation pattern: roll per-agent
  #      properties up into a per-node, per-tick report.
  set.seed(1)
  n <- 10000L; n_nodes <- 8L
  nid    <- sample.int(n_nodes, n, replace = TRUE) - 1L
  nodeid <- allocate_scalar("u16", n); nodeid$set(nid)
  census <- allocate_vector("i32", 5L, n_nodes)            # 5 ticks x 8 nodes
  bincount(nodeid, n_nodes, census, slot = 2L)             # write tick 2 (0-based)
  expect_equal(census$values()[3L, ], as.integer(tabulate(nid + 1L, n_nodes)))   # row 3 == slot 2
  expect_true(all(census$values()[1L, ] == 0L))            # other slots untouched
})
