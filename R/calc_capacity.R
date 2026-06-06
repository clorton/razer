# Estimate the agent capacity to preallocate for a growing population. Ported from
# laser.core's calc_capacity. For readers new to R: `inherits(x, "Column")` tests an
# object's S3 class; `^` is elementwise power on vectors/matrices; `colSums(m)` sums
# each column (here: over time, per node — numpy's `sum(axis=0)`); `stopifnot(cond)`
# errors if `cond` is not all TRUE; `.Machine$integer.max` is R's 32-bit signed int max.

#' Estimate the agent capacity needed for a growing population.
#'
#' Projects each node's population forward under its (possibly time-varying) crude birth
#' rate and returns a per-node capacity to preallocate, inflated by an optional safety
#' factor. The births are treated as geometric growth: a daily growth rate
#' `lambda = (1 + CBR/1000)^(1/365) - 1` is summed over all time steps and exponentiated
#' (`exp(sum lambda)`), a geometric-Brownian-motion-style estimate of expected growth.
#' The safety factor adds a multiple of the growth's standard-deviation-like term
#' `sqrt(exp(sum lambda)) - 1` as headroom for stochastic variation.
#'
#' @param birthrates A 2-D `nsteps x nnodes` numeric matrix of crude birth rates (births
#'   per 1,000 individuals per year), or a 2-D razer [Column] (e.g. from [values_map()])
#'   whose `$values()` is such a matrix. Each value must be in `[0, 100]`.
#' @param initial_pop A numeric vector of length `nnodes`: the initial population per
#'   node. Must be non-negative.
#' @param safety_factor Non-negative headroom multiplier in `[0, 6]` (default `1`). `0`
#'   gives the bare expected-growth estimate; larger values reserve more slack.
#' @return A numeric vector of length `nnodes` of estimated capacities (whole-valued
#'   doubles, which represent integers exactly up to `2^53`). A `warning` is issued if
#'   any estimate exceeds `.Machine$integer.max` (R's 32-bit signed integer max) — the
#'   largest count razer's allocators (`allocate_scalar`/`allocate_vector`, whose
#'   `count` is an `i32`) can accept.
#' @section Errors:
#' Stops if `birthrates` is not 2-D, if its node count does not match `length(initial_pop)`,
#' if any population is negative, if any birth rate is outside `[0, 100]`, or if
#' `safety_factor` is outside `[0, 6]`.
#' @examples
#' # two nodes, one year of a constant CBR of 40 per 1,000
#' br <- matrix(40, nrow = 365, ncol = 2)
#' calc_capacity(br, initial_pop = c(1e6, 5e5), safety_factor = 1)
#' @export
calc_capacity <- function(birthrates, initial_pop, safety_factor = 1) {
  # Accept a razer Column (a values_map grid) by copying out its matrix snapshot.
  if (inherits(birthrates, "Column")) birthrates <- birthrates$values()
  if (is.null(dim(birthrates)) || length(dim(birthrates)) != 2L)
    stop("`birthrates` must be a 2-D (nsteps x nnodes) matrix or a 2-D Column")

  initial_pop <- as.numeric(initial_pop)
  nnodes      <- ncol(birthrates)

  # ── validate (mirrors laser.core's assertions) ──────────────────────────────────
  if (length(initial_pop) != nnodes)
    stop(sprintf("number of nodes in `birthrates` (%d) and `initial_pop` length (%d) must match",
                 nnodes, length(initial_pop)))
  if (any(initial_pop < 0))
    stop("`initial_pop` values must be non-negative")
  if (any(birthrates < 0) || any(birthrates > 100))
    stop("all `birthrates` must be in [0, 100] (births per 1,000 per year)")
  if (!(length(safety_factor) == 1L && safety_factor >= 0 && safety_factor <= 6))
    stop(sprintf("`safety_factor` must be a single value in [0, 6], got %s", safety_factor))

  # ── projected growth ────────────────────────────────────────────────────────────
  # CBR (per 1,000 per year) -> per-individual daily growth rate. The annual growth
  # factor is (1 + CBR/1000); the daily factor is its 365th root, minus 1 for the rate.
  lamda <- (1 + birthrates / 1000)^(1 / 365) - 1
  # Expected growth factor per node: exp of the daily rates summed over all time steps
  # (geometric Brownian motion: E(P_t) = P_0 * exp(mu * t)).
  exp_mu_t <- exp(colSums(lamda))
  # Headroom: 1 + safety_factor * (sqrt(growth) - 1), a multiple of a sd-like term.
  safety_multiplier <- 1 + safety_factor * (sqrt(exp_mu_t) - 1)

  # Round to whole agents. Unlike laser.core (which returns a uint32 array and clamps to
  # uint32 max to avoid that dtype's modular wraparound), we return an R double, which
  # holds whole numbers exactly to 2^53 and never wraps — so no clamp is needed. We do
  # warn if any estimate exceeds R's 32-bit signed integer max, since that is the largest
  # count razer's allocators (allocate_scalar/allocate_vector, whose `count` is an i32)
  # can actually allocate.
  estimates <- round(initial_pop * safety_multiplier * exp_mu_t)
  over <- estimates > .Machine$integer.max
  if (any(over))
    warning(sprintf(
      "calc_capacity: %d node(s) exceed .Machine$integer.max (%d) and cannot be allocated as-is; largest is %.0f",
      sum(over), .Machine$integer.max, max(estimates)))
  estimates
}
