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
