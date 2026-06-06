# R-facing wrappers over the Rust bincount kernels. The four form a consistent family:
#   bincount()          count agents per bin
#   bincount_wt()       sum a weight  per bin
#   bincount_where()    count agents per bin, filtered by a predicate on a property
#   bincount_where_wt() sum a weight  per bin, filtered by a predicate on a property
# The `_wt` suffix adds weights; the `_where` suffix adds a predicate filter. The wrappers
# exist to give the `slot` argument a default of 0 (extendr's generated wrappers can't
# carry R-side defaults) and, for the `_where*` pair, to offer an allocate-and-return mode.
# For readers new to R: `function(args)` is a first-class closure bound with `<-`;
# `slot = 0L` is a default argument; `L` marks an integer literal; the last expression is
# the return value (the count kernels return `NULL` invisibly after mutating `counts`).

#' Count occurrences of each value, NumPy `bincount`-style, into a buffer.
#'
#' For each bin `b` in `0..nbins`, counts how many elements of `values` equal `b`
#' and writes the result into **slice `slot`** of `counts` — the `nbins` entries
#' `slot*slice_len .. slot*slice_len + nbins` — overwriting them and leaving the
#' rest untouched. For a scalar `counts` (shape `(1, n)`) the only slice is
#' `slot = 0`, the whole vector; for a 2-D report (e.g. `n_ticks x n_nodes` from
#' [allocate_vector()]) `slot` selects a tick's row. The tally is parallelized with
#' private per-thread histograms, so there are no write collisions.
#'
#' @param values An integer-typed [Column] of bin indices (`i8`..`u32`), each in
#'   `0..nbins`.
#' @param nbins  Number of bins; a non-negative integer no greater than `counts`'s
#'   slice length.
#' @param counts A numeric [Column] that receives the counts; modified in place.
#' @param slot   Which slice of `counts` to write; a non-negative integer less than
#'   `counts`'s slice count. Defaults to `0`.
#' @return `NULL`, invisibly; the result is written into `counts`.
#' @seealso [bincount_wt()], [bincount_where()], [bincount_where_wt()].
#' @examples
#' values <- allocate_scalar("u16", 6L)
#' values$set(c(0, 1, 1, 2, 2, 2))
#' counts <- allocate_scalar("i32", 3L)
#' bincount(values, 3L, counts)
#' counts$values()   # 1 2 3
#' @export
bincount <- function(values, nbins, counts, slot = 0L) {
  bincount_impl(values, nbins, counts, slot)
}

#' Weighted bincount: sum each element's weight into its bin.
#'
#' Like [bincount()], but accumulates `weights[i]` (rather than 1) into the bin
#' `values[i]`, à la `numpy.bincount(values, weights = ...)`. The per-bin sums are
#' written into **slice `slot`** of `counts`. `weights` must be the same length as
#' `values` and may be any numeric [Column] (signed, unsigned, or floating point).
#'
#' @param values  An integer-typed [Column] of bin indices (`i8`..`u32`).
#' @param weights A numeric [Column], the same length as `values`.
#' @param nbins   Number of bins; a non-negative integer no greater than `counts`'s
#'   slice length.
#' @param counts  A numeric [Column] that receives the weighted sums; modified in
#'   place.
#' @param slot    Which slice of `counts` to write; a non-negative integer less
#'   than `counts`'s slice count. Defaults to `0`.
#' @return `NULL`, invisibly; the result is written into `counts`.
#' @seealso [bincount()], [bincount_where()], [bincount_where_wt()].
#' @examples
#' values  <- allocate_scalar("u16", 5L); values$set(c(0, 0, 1, 2, 2))
#' weights <- allocate_scalar("f64", 5L); weights$set(c(1.5, 2.5, 4, 1, 3))
#' counts  <- allocate_scalar("f64", 3L)
#' bincount_wt(values, weights, 3L, counts)
#' counts$values()   # 4 4 4
#' @export
bincount_wt <- function(values, weights, nbins, counts, slot = 0L) {
  bincount_wt_impl(values, weights, nbins, counts, slot)
}

#' Count, per group, the agents whose property satisfies a comparison.
#'
#' A predicate-filtered, count-aware [bincount()]: for each group `g` in
#' `0..n_groups`, counts how many of the first `count` agents both have
#' `group[i] == g` AND satisfy `prop[i] <op> value`. This answers flexible
#' agent queries directly on the Columns — e.g. "exposed by node"
#' (`prop = state`, `op = "eq"`, `value = laser_states()[["E"]]`) or "under-fives
#' by node" (`prop = dob`, `op = "gt"`, `value = tick - 5 * 365`, since `dob` is the
#' negative age) — in one parallel pass with no copy of `prop` into R.
#'
#' Two output modes: leave `counts` `NULL` (the default) for an ad-hoc query and an
#' integer vector of length `n_groups` is allocated and returned; or pass a numeric
#' [Column] `counts` (e.g. a `n_ticks x n_nodes` report) and the totals are written
#' into its slice `slot` (and `NULL` returned invisibly), avoiding an allocation in a
#' per-tick model loop.
#'
#' @param group    An integer-typed [Column] of group indices (`i8`..`u32`), each in
#'   `0..n_groups` — typically `nodeid`.
#' @param n_groups Number of groups; a non-negative integer.
#' @param prop     A numeric [Column] holding the per-agent property to test
#'   (compared as a double).
#' @param op       Comparison string: one of `"eq"`, `"ne"`, `"lt"`, `"le"`, `"gt"`,
#'   `"ge"`.
#' @param value    The threshold the property is compared against.
#' @param count    How many leading agents to scan (the active count). Defaults to the
#'   full length of `group`.
#' @param counts   Optional numeric [Column] to receive the totals; when omitted an
#'   integer result vector is allocated and returned instead.
#' @param slot     Which slice of `counts` to write when `counts` is supplied; a
#'   non-negative integer less than `counts`'s slice count. Defaults to `0`.
#' @return When `counts` is `NULL`, an integer vector of per-group counts; otherwise
#'   `NULL` invisibly (the result is written into `counts`).
#' @seealso [bincount()], [bincount_wt()], [bincount_where_wt()].
#' @examples
#' states <- laser_states()
#' state  <- allocate_scalar("u8",  6L); state$set(c(0, 1, 1, 2, 1, 0))   # S E E I E S
#' nodeid <- allocate_scalar("u16", 6L); nodeid$set(c(0, 0, 1, 1, 1, 0))
#' bincount_where(nodeid, 2L, state, "eq", states[["E"]])   # exposed per node: 1 2
#' @export
bincount_where <- function(group, n_groups, prop, op, value,
                           count = NULL, counts = NULL, slot = 0L) {
  if (is.null(count)) count <- group$length()
  if (is.null(counts)) {
    # Ad-hoc mode: allocate a 1-D i32 result and hand back a plain integer vector.
    out <- allocate_scalar("i32", n_groups)
    bincount_where_impl(group, n_groups, prop, op, value, count, out, 0L)
    out$values()
  } else {
    # Report mode: write into the caller's buffer (slice `slot`); return NULL.
    bincount_where_impl(group, n_groups, prop, op, value, count, counts, slot)
  }
}

#' Sum a weight, per group, over the agents whose property satisfies a comparison.
#'
#' The weighted twin of [bincount_where()]: for each group `g` in `0..n_groups`, sums
#' `weights[i]` over the first `count` agents that both have `group[i] == g` AND satisfy
#' `prop[i] <op> value`. Use it for predicate-filtered weighted aggregates by node — e.g.
#' total infectiousness of symptomatic agents per node, or person-days under five — in one
#' parallel pass with no copy of `prop` or `weights` into R.
#'
#' Two output modes, as in [bincount_where()]: leave `counts` `NULL` (the default) and a
#' numeric vector of length `n_groups` is allocated and returned; or pass a numeric
#' [Column] `counts` and the per-group sums are written into its slice `slot`.
#'
#' @param group    An integer-typed [Column] of group indices (`i8`..`u32`) — e.g. `nodeid`.
#' @param n_groups Number of groups; a non-negative integer.
#' @param prop     A numeric [Column] holding the per-agent property to test.
#' @param op       Comparison string: one of `"eq"`, `"ne"`, `"lt"`, `"le"`, `"gt"`, `"ge"`.
#' @param value    The threshold the property is compared against.
#' @param weights  A numeric [Column] of per-agent weights to sum (any numeric type).
#' @param count    How many leading agents to scan (the active count). Defaults to the
#'   full length of `group`.
#' @param counts   Optional numeric [Column] to receive the sums; when omitted a numeric
#'   result vector is allocated and returned instead.
#' @param slot     Which slice of `counts` to write when `counts` is supplied. Defaults to `0`.
#' @return When `counts` is `NULL`, a numeric vector of per-group sums; otherwise `NULL`
#'   invisibly (the result is written into `counts`).
#' @seealso [bincount()], [bincount_wt()], [bincount_where()].
#' @examples
#' states  <- laser_states()
#' state   <- allocate_scalar("u8",  5L); state$set(c(2, 2, 0, 2, 1))     # I I S I E
#' nodeid  <- allocate_scalar("u16", 5L); nodeid$set(c(0, 0, 1, 1, 1))
#' shed    <- allocate_scalar("f64", 5L); shed$set(c(1.0, 0.5, 9, 2.0, 9))# infectiousness
#' # total shedding of infectious (state==I) agents per node:
#' bincount_where_wt(nodeid, 2L, state, "eq", states[["I"]], shed)        # 1.5  2.0
#' @export
bincount_where_wt <- function(group, n_groups, prop, op, value, weights,
                              count = NULL, counts = NULL, slot = 0L) {
  if (is.null(count)) count <- group$length()
  if (is.null(counts)) {
    out <- allocate_scalar("f64", n_groups)
    bincount_where_wt_impl(group, n_groups, prop, op, value, weights, count, out, 0L)
    out$values()
  } else {
    bincount_where_wt_impl(group, n_groups, prop, op, value, weights, count, counts, slot)
  }
}
