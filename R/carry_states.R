# Convenience wrapper that carries a set of per-node census counters forward one
# tick and (optionally) sums some of them into a running total — e.g. carry S, I, R
# forward and total them into N (the population) for the FOI denominator. For readers
# new to R: `lapply(xs, f)` maps `f` over a list; `Reduce(\`+\`, ...)` folds a list of
# numeric vectors with elementwise addition; `is.null` tests for the absent default.

#' Carry census counters forward, and optionally total some of them.
#'
#' For each Column in `carry`, calls [carry_forward()] (copying column `tick` onto
#' `tick + 1`). If `total` is supplied, then after carrying it sets `total`'s column
#' `tick + 1` to the elementwise sum of the `summands` Columns at `tick + 1` — for
#' example carrying `S`, `I`, `R` forward and totalling them into `N` so the current
#' per-node population is available to [calc_foi()] (and stays correct as births,
#' deaths, and imports change the compartments).
#'
#' @param carry     A list of 2-D census [Column]s to carry forward.
#' @param tick      0-based source tick; column `tick` is copied onto `tick + 1`.
#' @param total     Optional [Column] to receive the running total at column `tick + 1`.
#' @param summands  List of Columns to sum into `total` (defaults to `carry`).
#' @return `NULL`, invisibly; the Columns are modified in place.
#' @examples
#' \dontrun{
#' # carry S, I, R forward and keep N = S + I + R up to date:
#' carry_forward_states(list(nodes$S, nodes$I, nodes$R), tick, total = nodes$N)
#' }
#' @export
carry_forward_states <- function(carry, tick, total = NULL, summands = carry) {
  for (col in carry) carry_forward(col, tick)
  if (!is.null(total)) {
    next_col <- tick + 1L
    parts <- lapply(summands, function(col) col$col(next_col))
    total$set_col(next_col, Reduce(`+`, parts))
  }
  invisible(NULL)
}
