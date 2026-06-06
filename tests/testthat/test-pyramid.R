# Tests for the R-side pyramid helpers: load_pyramid_csv() parses a population-pyramid
# CSV into a start/end/M/F matrix, and sample_pyramid_ages() draws per-agent ages (in
# days) from it via the alias sampler. Written given-when-then. `withr`-free: we write a
# tempfile() and clean up with on.exit.

valid_csv <- c(
  "Age,M,F",
  "0-4,100,98",
  "5-9,90,92",
  "10+,5,6"
)

test_that("load_pyramid_csv parses bands into a start/end/M/F matrix", {
  # Given a well-formed pyramid CSV with two ranges and an open-ended final band
  # When it is loaded
  # Then the result is a 3x4 integer matrix with the expected columns, and the final
  #      "10+" band is stored as a single-year bucket (end == start == 10). Failure
  #      would mean the parser mishandles the range split or the open-ended last band.
  path <- tempfile(fileext = ".csv"); on.exit(unlink(path))
  writeLines(valid_csv, path)

  pyr <- load_pyramid_csv(path)

  expect_equal(colnames(pyr), c("start", "end", "M", "F"))
  expect_equal(unname(pyr[1L, ]), c(0L, 4L, 100L, 98L))
  expect_equal(unname(pyr[2L, ]), c(5L, 9L, 90L, 92L))
  expect_equal(unname(pyr[3L, ]), c(10L, 10L, 5L, 6L))   # "10+" -> single-year bucket
})

test_that("load_pyramid_csv rejects malformed files", {
  # Given files with a bad header, a malformed data line, and non-ascending ages
  # When each is loaded
  # Then load_pyramid_csv errors with a message identifying the problem.
  # Failure would mean bad demographic input is silently accepted.
  bad_header <- tempfile(fileext = ".csv"); on.exit(unlink(bad_header), add = TRUE)
  writeLines(c("Age,Male,Female", "0-4,1,1", "5+,1,1"), bad_header)
  expect_error(load_pyramid_csv(bad_header), "Age,M,F")

  bad_line <- tempfile(fileext = ".csv"); on.exit(unlink(bad_line), add = TRUE)
  writeLines(c("Age,M,F", "0to4,1,1", "5+,1,1"), bad_line)
  expect_error(load_pyramid_csv(bad_line), "low-high")

  unsorted <- tempfile(fileext = ".csv"); on.exit(unlink(unsorted), add = TRUE)
  writeLines(c("Age,M,F", "5-9,1,1", "0-4,1,1", "10+,1,1"), unsorted)
  expect_error(load_pyramid_csv(unsorted), "ascending")
})

test_that("sample_pyramid_ages draws ages within the pyramid's day range", {
  # Given a small pyramid (bands [0,4], [5,9], [10,10] years) where the youngest band
  #       dominates the population
  # When 50,000 ages (in days) are drawn
  # Then every age is in [0, 11*365) days, and the youngest band [0,5) years is the
  #      most populous in the sample (it has the largest count). Failure would mean the
  #      band-to-age-range mapping or the population weighting is wrong.
  set.seed(1L)
  pyr <- rbind(c(0, 4, 1000, 1000), c(5, 9, 100, 100), c(10, 10, 10, 10))
  colnames(pyr) <- c("start", "end", "M", "F")

  ages <- sample_pyramid_ages(pyr, 50000L)

  expect_length(ages, 50000L)
  expect_true(all(ages >= 0L & ages < 11L * 365L))
  # most agents should be under 5 years (the heaviest band)
  expect_gt(mean(ages < 5L * 365L), 0.8)
})
