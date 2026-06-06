# R-facing wrappers over the Rust `bincount_impl` / `bincountw_impl` kernels, added
# solely to give the `slot` argument a default of 0. extendr's generated wrappers
# (in R/extendr-wrappers.R) can't carry R-side default arguments, so these thin
# wrappers supply one. For readers new to R: `function(args)` is a first-class
# closure value bound with `<-`; `slot = 0L` is a default argument (used when the
# caller omits it); `L` marks an integer literal; the last expression is the return
# value (here the kernels return `NULL` invisibly after mutating `counts`).

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
#' @examples
#' values  <- allocate_scalar("u16", 5L); values$set(c(0, 0, 1, 2, 2))
#' weights <- allocate_scalar("f64", 5L); weights$set(c(1.5, 2.5, 4, 1, 3))
#' counts  <- allocate_scalar("f64", 3L)
#' bincountw(values, weights, 3L, counts)
#' counts$values()   # 4 4 4
#' @export
bincountw <- function(values, weights, nbins, counts, slot = 0L) {
  bincountw_impl(values, weights, nbins, counts, slot)
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
#' @examples
#' states <- laser_states()
#' state  <- allocate_scalar("u8",  6L); state$set(c(0, 1, 1, 2, 1, 0))   # S E E I E S
#' nodeid <- allocate_scalar("u16", 6L); nodeid$set(c(0, 0, 1, 1, 1, 0))
#' count_by_where(nodeid, 2L, state, "eq", states[["E"]])   # exposed per node: 1 2
#' @export
count_by_where <- function(group, n_groups, prop, op, value,
                           count = NULL, counts = NULL, slot = 0L) {
  if (is.null(count)) count <- group$length()
  if (is.null(counts)) {
    # Ad-hoc mode: allocate a 1-D i32 result and hand back a plain integer vector.
    out <- allocate_scalar("i32", n_groups)
    count_by_where_impl(group, n_groups, prop, op, value, count, out, 0L)
    out$values()
  } else {
    # Report mode: write into the caller's buffer (slice `slot`); return NULL.
    count_by_where_impl(group, n_groups, prop, op, value, count, counts, slot)
  }
}
