# squash() — compact a people environment, reclaiming the slots of agents the mask
# excludes (by default, the deceased). The Column model keeps each agent property in a
# separate capacity-sized Column; squash() applies ONE keep-mask to every per-agent
# Column so they stay aligned, then updates the active count. The per-node census (in the
# `nodes` environment) is aggregate and unaffected — only the agent arrays are compacted.
# For readers new to R: `ls(env)` lists an environment's names; `get(nm, envir = env)`
# fetches one; `inherits(x, "Column")` is an isinstance check.

#' Compact a people environment, dropping excluded agents and reclaiming their slots.
#'
#' Applies a logical `keep` mask (length `people$count`) to every per-agent [Column] in
#' the `people` environment — shifting the survivors to the front of each array, in order
#' — and sets `people$count` to the number kept. Reuse frees the slots of agents that have
#' left the simulation (e.g. the deceased) so the per-tick kernels stop iterating over
#' them and `births` can refill the slots. All Columns are compacted by the SAME mask, so
#' they remain row-aligned.
#'
#' @param people An environment whose per-agent properties are scalar [Column]s (e.g.
#'   `state`, `timer`, `nodeid`, `dob`, `dod`) plus an integer `count`.
#' @param keep Optional logical vector of length `people$count`; `TRUE` keeps the agent.
#'   Defaults to "still alive" — every agent whose `state` is not `D` (stored as 255 in
#'   the u8 state Column).
#' @return The new active count (invisibly); `people$count` is updated in place.
#' @examples
#' \dontrun{
#' # Periodically reclaim dead agents during a long run with mortality:
#' people$count <- squash(people)
#' }
#' @export
squash <- function(people, keep = NULL) {
  n <- people$count
  if (is.null(keep)) {
    # D (-1) is stored as 255 in the u8 `state` Column; everything else is alive.
    keep <- people$state$values()[seq_len(n)] != 255L
  } else {
    keep <- as.logical(keep)
    if (length(keep) != n)
      stop(sprintf("`keep` must have length people$count (%d), got %d", n, length(keep)))
  }
  n_keep <- n
  for (nm in ls(people)) {
    col <- get(nm, envir = people)
    if (inherits(col, "Column")) n_keep <- col$squash(keep)   # same mask -> stays aligned
  }
  people$count <- n_keep
  invisible(n_keep)
}
