# =============================================================================
# models.R — high-level model *runners* that assemble the per-tick step kernels
# (epidemic.rs) and their parameters into a single call.
#
# `run_sir()` is the reference template; `run_seir()`, `run_si()`, `run_sis()`, …
# will follow exactly the same shape. Each runner:
#   1. takes a data.frame describing the geography (one row per node, an integer
#      `population` column) — R's analogue of a pandas DataFrame;
#   2. builds an agent `LaserFrame` ("people") sized to the total population and
#      assigns every agent to its node;
#   3. seeds initial infections, then runs the tick loop, calling the kernels in
#      downstream-first order (see CLAUDE.md), recording per-node compartment
#      trajectories and inter-compartment flows; and
#   4. returns the people frame, with the recorded time series attached as
#      attributes.
#
# Orientation for readers coming from C / C++ / C# / Python (see laser_frame.R
# for the longer primer):
#   * `<-` is assignment; `function(args) body` is a first-class closure; the
#     last expression in a function is its return value (no `return` needed).
#   * `pkg::fn` is a namespaced call; `x[["k"]]` extracts one named element.
#   * R is 1-based: `seq_len(k)` is the vector `1, 2, …, k`.
#   * `rep(x, times = n)` repeats elements; `cumsum` is a running total.
#   * `tabulate(v, nbins)` is a fixed-width histogram — it counts how many of
#     `1..nbins` occur in `v` (values outside the range are ignored), always
#     returning a length-`nbins` vector (so empty nodes report 0, not a gap).
#   * Arithmetic/comparison are vectorized: `state == code` yields a logical
#     vector, and `a - b` subtracts element-wise — no explicit loop required.
#   * `attr(x, "k") <- v` attaches named metadata to an object without changing
#     its class; `inherits(x, "Cls")` is R's `isinstance` check.
# =============================================================================

# Coerce a user-supplied duration argument into a `Distribution`. A bare number
# is promoted to a constant (fixed-duration) distribution, so a caller can pass
# either `7` or `dist_gamma(2, 3.5)` interchangeably. Internal helper (no @export).
as_distribution <- function(x, what = "argument") {
  if (inherits(x, "Distribution")) return(x)
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) return(dist_constant(x))
  stop(sprintf(
    "`%s` must be a Distribution (e.g. dist_gamma(...)) or a single finite number, not %s",
    what, paste(class(x), collapse = "/")
  ))
}

# Round duration draws to whole ticks, clamped to a minimum of 1 — the R-side
# mirror of epidemic.rs::duration_ticks. Used only to set the seeded agents'
# initial infectious timers; the kernels apply the same rule internally for
# every later transition. `pmax` is the element-wise (parallel) maximum.
.duration_ticks <- function(x) pmax(as.integer(round(x)), 1L)

# Count active agents per node that are currently in compartment `code`.
# `state == code` is a vectorized logical mask; `node_ids[mask]` keeps the
# 0-based node id of each matching agent; `+ 1L` shifts to 1-based so tabulate's
# `nbins = n_nodes` pins the result length to the node count.
count_by_node <- function(state, node_ids, code, n_nodes) {
  tabulate(node_ids[state == code] + 1L, nbins = n_nodes)
}

#' Run an agent-based SIR model.
#'
#' Builds a population of individual agents from a per-node scenario table, seeds
#' initial infections, and advances a stochastic SIR model for `nticks` ticks.
#' Each tick applies the kernels downstream-first — recovery (I→R) before
#' transmission (S→I) — so a newly infectious agent is not decremented in the
#' same tick it is infected (see the modeling-convention notes in `CLAUDE.md`).
#'
#' The returned people `LaserFrame` holds the *final* per-agent state. The
#' recorded time series are attached as attributes:
#'
#' * `attr(model, "report")` — a node-level `LaserFrame` (one row per node) whose
#'   vector properties are the recorded matrices, each `n_nodes` rows by time:
#'   - `S`, `I`, `R`: compartment counts at every tick boundary (`nticks + 1`
#'     columns; column 1 is the initial state, column `t + 1` is the state after
#'     tick `t`).
#'   - `incidence` (S→I) and `recovery` (I→R): per-tick flow counts (`nticks`
#'     columns; column `t` covers tick `t`).
#'
#'   Pull a matrix with `attr(model, "report")$I` (counts) or transpose with
#'   `t(...)` to get time-down-the-rows.
#' * `attr(model, "runtime")` — wall-clock seconds spent in the tick loop.
#' * `attr(model, "nodes")` — the node frame passed to the transmission kernel.
#' * `attr(model, "nticks")`, `attr(model, "model")`, `attr(model, "parameters")`.
#'
#' @param nodes A data.frame with one row per geographic node and an integer
#'   `population` column giving each node's initial population.
#' @param beta Transmission rate (force of infection per infectious contact per
#'   tick).
#' @param infectious_period Infectious period in ticks: either a `Distribution`
#'   (e.g. `dist_gamma(...)`) or a single number (a fixed period).
#' @param nticks Number of ticks to simulate (a positive integer).
#' @param seed Initial number of infectious agents per node: either one number
#'   (applied to every node) or a vector of length `nrow(nodes)`. Defaults to 0
#'   (no infections — supply a positive seed to start an epidemic).
#' @param progress If `TRUE`, draw a text progress bar during the run.
#'
#' @return The people `LaserFrame` (final per-agent state) with the recorded
#'   trajectories and flows attached as attributes (see Details).
#'
#' @examples
#' nodes <- data.frame(population = c(10000, 5000))
#' model <- run_sir(nodes, beta = 0.3, infectious_period = 7, nticks = 160,
#'                  seed = c(10, 0))
#' report <- attr(model, "report")
#' tail(t(report$I))        # infectious count over time, per node
#' colSums(report$incidence)  # total new infections per node
#' attr(model, "runtime")
#' @export
run_sir <- function(nodes, beta, infectious_period, nticks,
                    seed = 0L, progress = FALSE) {
  # ── validate the scenario table ─────────────────────────────────────────────
  if (!is.data.frame(nodes))
    stop("`nodes` must be a data.frame with one row per node")
  if (!"population" %in% names(nodes))
    stop("`nodes` must have an integer `population` column")
  pops <- as.integer(nodes$population)
  if (anyNA(pops) || any(pops < 0L))
    stop("`population` must be non-negative integers")
  n_nodes <- nrow(nodes)
  n <- sum(pops)
  if (n <= 0L) stop("total population must be positive")

  nticks <- as.integer(nticks)
  if (length(nticks) != 1L || is.na(nticks) || nticks < 1L)
    stop("`nticks` must be a positive integer")

  inf_dist <- as_distribution(infectious_period, "infectious_period")

  # Per-node seed counts: a single value is recycled to every node.
  seed <- as.integer(seed)
  if (length(seed) == 1L) seed <- rep(seed, n_nodes)
  if (length(seed) != n_nodes)
    stop("`seed` must have length 1 or nrow(nodes)")
  if (anyNA(seed) || any(seed < 0L) || any(seed > pops))
    stop("`seed` must be between 0 and each node's population")

  # State codes as a named integer vector: states[["S"]] == 0L, etc.
  states <- laser_states()
  S <- states[["S"]]; I <- states[["I"]]; R <- states[["R"]]

  # ── build the agent (people) frame ──────────────────────────────────────────
  people <- LaserFrame$new(n, n)
  people$add_scalar_property("state", "integer", S)    # everyone starts susceptible
  people$add_scalar_property("node",  "integer", 0L)
  people$add_scalar_property("timer", "integer", 0L)

  # Assign agents to nodes: node 0 gets the first pops[1] agents, node 1 the next
  # pops[2], and so on. `seq_len(n_nodes) - 1L` are the 0-based node codes the
  # kernels expect; `rep(..., times = pops)` expands them to one entry per agent.
  node_ids <- rep(seq_len(n_nodes) - 1L, times = pops)
  people$node <- node_ids                              # `$<-` writes the whole column

  # ── seed initial infections per node ────────────────────────────────────────
  # offsets[k] is the 0-based index of node k's first agent, so node k's agents
  # occupy positions (offsets[k] + 1) .. offsets[k + 1] within the frame.
  offsets <- c(0L, cumsum(pops))
  state0 <- people$state                               # length-n vector, all S
  timer0 <- people$timer                               # length-n vector, all 0
  for (k in seq_len(n_nodes)) {
    n_seed <- seed[k]
    if (n_seed > 0L) {
      idx <- (offsets[k] + 1L):(offsets[k] + n_seed)
      state0[idx] <- I
      # Seeded agents get an infectious timer drawn from the same distribution
      # the kernel uses for newly infected agents (rounded, clamped to >= 1).
      timer0[idx] <- .duration_ticks(inf_dist$sample_n(n_seed))
    }
  }
  people$state <- state0
  people$timer <- timer0

  # Node frame consumed by step_transmission_si: it needs the per-node population
  # `N` and an `I` column the kernel overwrites with the current infectious count.
  nd <- LaserFrame$new(n_nodes, n_nodes)
  nd$add_scalar_property("N", "integer", 0L)
  nd$N <- pops
  nd$add_scalar_property("I", "integer", 0L)

  # SIR: R is absorbing, so the recovered timer is never read; dist_constant(0)
  # is the conventional "don't care" duration for step_infectious_ir.
  imm_dist <- dist_constant(0)

  # ── allocate the report (one row per node, time across the columns) ──────────
  report <- LaserFrame$new(n_nodes, n_nodes)
  # Compartment trajectories: nticks + 1 columns (column 1 = initial state).
  report$add_vector_property("S", nticks + 1L, "integer", 0L)
  report$add_vector_property("I", nticks + 1L, "integer", 0L)
  report$add_vector_property("R", nticks + 1L, "integer", 0L)
  # Flows: one count per tick (nticks columns).
  report$add_vector_property("incidence", nticks, "integer", 0L)  # S -> I
  report$add_vector_property("recovery",  nticks, "integer", 0L)  # I -> R

  # Snapshot per-node compartment counts from the current agent states. One pull
  # of the state column (a copy from Rust), then three histograms over it. This
  # closure captures `people`, `node_ids`, `n_nodes`, and the S/I/R codes.
  tally <- function() {
    st <- people$state
    list(S = count_by_node(st, node_ids, S, n_nodes),
         I = count_by_node(st, node_ids, I, n_nodes),
         R = count_by_node(st, node_ids, R, n_nodes))
  }

  # Record the initial state into column 1 of each trajectory.
  cur <- tally()
  report$set_col("S", 1L, cur$S)
  report$set_col("I", 1L, cur$I)
  report$set_col("R", 1L, cur$R)

  # `txtProgressBar` returns a handle updated with `setTxtProgressBar`; NULL means
  # "no bar". `isTRUE` guards against a non-logical `progress` argument.
  pb <- if (isTRUE(progress))
    utils::txtProgressBar(min = 0L, max = nticks, style = 3L) else NULL

  # ── tick loop ───────────────────────────────────────────────────────────────
  # Flows are measured as the drop in a *source* compartment around its step:
  # each kernel drains exactly one compartment with no inflow during that step,
  # so the per-node count before minus after is exactly the number that moved.
  # `system.time({ ... })` runs the block and returns a timing vector; the
  # assignments inside still take effect in this scope. `[["elapsed"]]` is the
  # wall-clock field.
  elapsed <- system.time({
    for (tick in seq_len(nticks)) {
      # Downstream-first: drain I -> R before S -> I.
      step_infectious_ir(people, imm_dist = imm_dist)
      mid <- tally()                       # after recovery, before transmission
      recovery <- cur$I - mid$I            # agents that left I this tick, per node

      step_transmission_si(people, nd, beta = beta, inf_dist = inf_dist)
      new <- tally()                       # after transmission
      incidence <- cur$S - new$S           # agents that left S this tick, per node

      report$set_col("recovery",  tick, recovery)
      report$set_col("incidence", tick, incidence)
      report$set_col("S", tick + 1L, new$S)
      report$set_col("I", tick + 1L, new$I)
      report$set_col("R", tick + 1L, new$R)

      cur <- new                           # carry forward for the next tick's deltas
      if (!is.null(pb)) utils::setTxtProgressBar(pb, tick)
    }
  })[["elapsed"]]
  if (!is.null(pb)) close(pb)

  # ── attach the recorded series + run metadata to the returned people frame ───
  attr(people, "report")     <- report
  attr(people, "nodes")      <- nd
  attr(people, "node_ids")   <- node_ids
  attr(people, "model")      <- "SIR"
  attr(people, "nticks")     <- nticks
  attr(people, "runtime")    <- elapsed
  attr(people, "parameters") <- list(beta = beta, infectious_period = inf_dist)
  people
}
