# Population-pyramid helpers: read a pyramid CSV, and turn it into realistic per-agent
# ages using the Vose alias sampler (`aliased_distribution`, from the Rust side). For
# readers new to R: `readLines` returns a character vector (one element per line);
# `trimws` strips surrounding whitespace; `grepl(pattern, x)` returns a logical vector
# of regex matches; `strsplit(x, sep)` returns a LIST of split pieces; `[[1]]` indexes
# the first list element; `do.call(rbind, list_of_rows)` stacks rows into a matrix.

#' Load a population-pyramid CSV into a numeric matrix.
#'
#' The file is expected to have the schema used by laser.core pyramids:
#'
#' - a header line exactly `"Age,M,F"`;
#' - then one line per age band `"low-high,males,females"` (all non-negative integers);
#' - a final open-ended band `"max+,males,females"`.
#'
#' The returned matrix has one row per band and the integer columns `start`, `end`,
#' `M`, `F`. The final `"max+"` band is stored as a single-year bucket (`end == start`),
#' matching laser.core.
#'
#' @param file Path to the CSV file.
#' @return An integer matrix with columns `start`, `end`, `M`, `F` (one row per band).
#' @section Errors:
#' Stops if the header is not `"Age,M,F"`, if any data line is malformed, or if the
#' start/end ages are not strictly ascending.
#' @examples
#' \dontrun{
#' pyramid <- load_pyramid_csv("USA_pyramid_2020.csv")
#' }
#' @export
load_pyramid_csv <- function(file) {
  lines <- trimws(readLines(file))               # one trimmed string per line
  if (length(lines) < 2L)
    stop("pyramid file must have a header and at least one data line")
  if (lines[1L] != "Age,M,F")
    stop("header line is not 'Age,M,F'")

  body <- lines[-1L]                             # drop the header
  n    <- length(body)
  # All but the last band match "low-high,males,females"; the last matches "max+,…".
  mid_ok  <- grepl("^[0-9]+-[0-9]+,[0-9]+,[0-9]+$", body[-n])
  last_ok <- grepl("^[0-9]+\\+,[0-9]+,[0-9]+$", body[n])
  if (!all(mid_ok))
    stop("data lines must look like 'low-high,males,females'")
  if (!last_ok)
    stop("last data line must look like 'max+,males,females'")

  # Parse each line into c(start, end, M, F). `as.integer` on the split pieces; the
  # final band's "max+" has no "-high", so its end is set equal to its start.
  parse_row <- function(line, last) {
    parts <- strsplit(line, ",")[[1L]]           # e.g. c("0-4", "100", "98")
    counts <- as.integer(parts[2:3])             # males, females
    if (last) {
      start <- as.integer(sub("\\+$", "", parts[1L]))
      c(start, start, counts)
    } else {
      ages <- as.integer(strsplit(parts[1L], "-")[[1L]])  # c(low, high)
      c(ages, counts)
    }
  }
  rows <- lapply(seq_len(n), function(i) parse_row(body[i], i == n))
  mat  <- do.call(rbind, rows)                   # n x 4 integer matrix
  colnames(mat) <- c("start", "end", "M", "F")

  # Validity: starting and ending ages strictly ascending (as in laser.core).
  if (!all(diff(mat[, "start"]) > 0L))
    stop("starting ages are not in ascending order")
  if (!all(diff(mat[, "end"]) > 0L))
    stop("ending ages are not in ascending order")

  mat
}

#' Sample realistic per-agent ages (in days) from a population pyramid.
#'
#' Builds an [aliased_distribution()] over the per-band populations (`M + F`), draws a
#' band for each agent in proportion to its population, then a uniform day within that
#' band's year range `[start, end + 1)` — so an agent in the `0-4` band gets an age
#' uniformly in `[0, 5)` years.
#'
#' @param pyramid An integer matrix with columns `start`, `end`, `M`, `F`, as returned
#'   by [load_pyramid_csv()].
#' @param n Number of agent ages to draw.
#' @return An integer vector of length `n` of ages in whole days.
#' @details Band selection uses the package's internal (thread-local, not R-seedable)
#'   RNG via the alias sampler; the within-band day uses R's RNG (`set.seed`-able).
#' @examples
#' pyramid <- rbind(c(0, 4, 100, 98), c(5, 9, 90, 92), c(10, 10, 5, 6))
#' colnames(pyramid) <- c("start", "end", "M", "F")
#' ages_days <- sample_pyramid_ages(pyramid, 1000L)
#' @export
sample_pyramid_ages <- function(pyramid, n) {
  n      <- as.integer(n)
  counts <- as.numeric(pyramid[, "M"] + pyramid[, "F"])
  dist   <- aliased_distribution(counts)
  bins   <- dist$sample_n(n) + 1L                # 0-based bins -> 1-based matrix rows
  lo     <- pyramid[bins, "start"] * 365L        # band lower bound, in days
  hi     <- (pyramid[bins, "end"] + 1L) * 365L   # band upper bound (exclusive), in days
  as.integer(lo + floor(stats::runif(n) * (hi - lo)))
}
