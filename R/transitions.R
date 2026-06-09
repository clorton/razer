# Apply a per-node transition count returned by a kernel to the node census. The
# Column kernels (transmission, step_si/step_sir/step_sirs, mortality, births) mutate the
# per-agent arrays and RETURN per-node counts; the model uses move_count() to push those
# counts into whichever census buffers it maintains. For readers new to R: `is.null(x)`
# tests for the absent default; `col$col(t)` / `col$set_col(t, v)` read / write one column
# (one tick) of a 2-D Column.

#' Apply a per-node transition count to the census at `tick + 1`.
#'
#' Subtracts `counts` from the `from` state and adds it to the `to` state at
#' census column `tick + 1` (the working column the model has already carried forward).
#' Either side may be `NULL` to skip it — e.g. a death is a one-sided decrement
#' (`to = NULL`) and a birth a one-sided increment (`from = NULL`).
#'
#' @param from A 2-D census [Column] to decrement, or `NULL`.
#' @param to   A 2-D census [Column] to increment, or `NULL`.
#' @param counts Integer vector of per-node counts (length `n_nodes`).
#' @param tick 0-based source tick; the delta is applied at column `tick + 1`.
#' @return `NULL`, invisibly; the Columns are modified in place.
#' @examples
#' \dontrun{
#' inf <- transmission(state, timer, nodeid, count, nodes$foi, t, states[["E"]], inc_dur)
#' move_count(nodes$S, nodes$E, inf, t)   # S -> E
#' }
#' @export
move_count <- function(from, to, counts, tick) {
  t1 <- tick + 1L
  if (!is.null(from)) from$set_col(t1, from$col(t1) - counts)
  if (!is.null(to))   to$set_col(t1,   to$col(t1)   + counts)
  invisible(NULL)
}
