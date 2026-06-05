# Utility for building a "values map" — a time x node parameter grid as a Column,
# used for model inputs that may vary over time and/or space (e.g. the transmission
# coefficient beta and a seasonal modifier passed to calc_foi()). The name matches
# LASER's ValuesMap. For readers new to R: `<-` is assignment; `is.matrix`/`length`/
# `dim` inspect the argument's shape; `rep(x, times)` repeats the whole vector,
# `rep(x, each)` repeats each element; `t()` transposes; `stop()` raises an error.

#' Build a values map: broadcast a value into a `n_ticks x n_nodes` grid [Column].
#'
#' Expands a flexible per-time and/or per-node value into a full `n_ticks x n_nodes`
#' f64 [Column] (slice-major: each tick's per-node row contiguous), suitable for
#' passing to [calc_foi()] as `beta` or `seasonality`. This is razer's equivalent of
#' LASER's *ValuesMap*. The shape of `value` selects how it is broadcast:
#'
#' * **scalar** (length 1) — constant over time and space.
#' * **length `n_nodes`** — varies by node, constant over time (per-node).
#' * **length `n_ticks`** — varies by time, constant over space (per-tick).
#' * **`n_ticks x n_nodes` matrix** — varies by both; used as-is.
#'
#' When `n_ticks == n_nodes` a bare vector is ambiguous; it is treated as per-node.
#' Pass an explicit matrix to vary by both in that case.
#'
#' @param value   A scalar, a length-`n_nodes` or length-`n_ticks` numeric vector,
#'   or a `n_ticks x n_nodes` numeric matrix.
#' @param n_ticks Number of time slices (rows of the grid).
#' @param n_nodes Number of nodes (columns of the grid).
#' @return An f64 [Column] of shape `n_ticks x n_nodes`.
#' @examples
#' g <- values_map(0.5, 10L, 3L)        # constant 0.5 everywhere
#' dim(g$values())                       # 10 3
#' values_map(c(1, 2, 3), 10L, 3L)       # per-node (length n_nodes)
#' values_map(seq_len(10L), 10L, 3L)     # per-tick (length n_ticks)
#' @export
values_map <- function(value, n_ticks, n_nodes) {
  n_ticks <- as.integer(n_ticks)
  n_nodes <- as.integer(n_nodes)
  # Build the buffer in SLICE-MAJOR order: element [tick, node] at tick*n_nodes + node.
  vec <-
    if (is.matrix(value)) {
      if (!all(dim(value) == c(n_ticks, n_nodes)))
        stop(sprintf("matrix `value` must be %d x %d, got %d x %d",
                     n_ticks, n_nodes, nrow(value), ncol(value)))
      as.numeric(t(value))                       # transpose -> tick-major
    } else if (length(value) == 1L) {
      rep(as.numeric(value), n_ticks * n_nodes)  # scalar
    } else if (length(value) == n_nodes) {
      rep(as.numeric(value), times = n_ticks)    # per-node, repeated each tick
    } else if (length(value) == n_ticks) {
      rep(as.numeric(value), each = n_nodes)     # per-tick, repeated across nodes
    } else {
      stop(sprintf(paste0("`value` must be a scalar, a length-%d (per-node) or ",
                          "length-%d (per-tick) vector, or a %d x %d matrix; got length %d"),
                   n_nodes, n_ticks, n_ticks, n_nodes, length(value)))
    }
  grid <- allocate_vector("f64", n_ticks, n_nodes)
  grid$set(vec)
  grid
}
