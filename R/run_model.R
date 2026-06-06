# run_model() — a high-level runner for the closed-population compartmental menagerie
# (SI/SEI/SIS/SEIS/SIR/SEIR/SIRS/SEIRS) on the Column kernels. It owns the per-tick
# bookkeeping — the carry-forward, the `calc_foi` placement that yields R0 = beta*D, and
# the move_count census deltas — so a user gets a correct model from one call instead of
# ~100 lines of hand-wiring (and there's no way to desync the census from the agents).
#
# For readers new to R: `switch`/`if` choose per-model behaviour; closures (`function(t)`)
# capture the people/nodes/distribution objects; `inherits(x,"Distribution")` is an
# isinstance check; `grepl("E", model)` tests for a substring.

# Coerce a duration argument to a Distribution (a bare number becomes a constant).
.as_dist <- function(x, what) {
  if (inherits(x, "Distribution")) return(x)
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) return(dist_constant(x))
  stop(sprintf("`%s` must be a Distribution or a single finite number", what))
}

# Round duration draws to whole u16 ticks (the kernels' minimum is 1).
.timer0 <- function(x) pmax(1L, pmin(65535L, as.integer(round(x))))

#' Run a closed-population compartmental model.
#'
#' Builds an agent population from a per-node scenario, seeds it, and advances one of the
#' eight SI / SEI / SIS / SEIS / SIR / SEIR / SIRS / SEIRS models for `nticks` ticks using
#' the Column kernels in the correct order (so the realized basic reproduction number is
#' the full `R0 = beta * D`). The model is `transmission` (or `transmission_si` for SI)
#' plus one step kernel; `run_model` applies every census delta for you.
#'
#' @param scenario A data.frame with one row per node and an integer `population` column;
#'   optional integer `I` / `R` columns seed initial infectious / recovered per node.
#' @param model One of `"SI"`, `"SEI"`, `"SIS"`, `"SEIS"`, `"SIR"`, `"SEIR"`, `"SIRS"`,
#'   `"SEIRS"` (case-insensitive).
#' @param nticks Number of recorded daily states (dynamics run `nticks - 1` times).
#' @param r0 Basic reproduction number; `beta = r0 / mean(infectious_period)`.
#' @param infectious_period Infectious period in ticks (a Distribution or a number).
#' @param incubation_period Latent/exposed period (required for the `SE*` models).
#' @param immunity_period Waning-immunity period (required for `SIRS` / `SEIRS`).
#' @param network Optional `n_nodes x n_nodes` spatial-coupling matrix (default: an
#'   all-zero matrix — independent nodes).
#' @param seasonality Transmission modifier in any [values_map()]-broadcastable form
#'   (default `1`).
#' @param seed Optional integer seed; if supplied, [set_seed()] is called for a
#'   reproducible run.
#' @param progress If `TRUE`, draw a text progress bar.
#' @return Invisibly, a list with the `people` and `nodes` environments. `nodes` holds the
#'   per-tick census Columns the model uses (`S`, and `E`/`I`/`R` as applicable), `N`, the
#'   `incidence` flow (new infections), and (for models with an I-exit) a `recoveries`
#'   flow — each an `n_nodes`-wide, time-major Column.
#' @examples
#' \dontrun{
#' scenario <- data.frame(population = 1e5L, I = 100L)
#' res <- run_model(scenario, "SEIR", nticks = 200L, r0 = 2.5,
#'                  infectious_period = 8, incubation_period = 5, seed = 1)
#' tail(res$nodes$I$values())
#' }
#' @export
run_model <- function(scenario, model, nticks, r0, infectious_period,
                      incubation_period = NULL, immunity_period = NULL,
                      network = NULL, seasonality = 1, seed = NULL, progress = FALSE) {
  model <- toupper(model)
  valid <- c("SI", "SEI", "SIS", "SEIS", "SIR", "SEIR", "SIRS", "SEIRS")
  if (!model %in% valid)
    stop(sprintf("`model` must be one of %s", paste(valid, collapse = ", ")))
  if (!is.data.frame(scenario) || !"population" %in% names(scenario))
    stop("`scenario` must be a data.frame with an integer `population` column")
  nticks <- as.integer(nticks)
  if (length(nticks) != 1L || is.na(nticks) || nticks < 2L)
    stop("`nticks` must be an integer >= 2")
  if (!is.null(seed)) set_seed(seed)

  states  <- laser_states()
  pops    <- as.integer(scenario$population)
  n       <- sum(pops)
  n_nodes <- nrow(scenario)

  has_E   <- grepl("E", model)
  has_R   <- grepl("R", model)
  waning  <- model %in% c("SIRS", "SEIRS")             # I->R sets an immunity timer, R->S
  has_step_clearance <- !model %in% c("SI", "SEI")     # I leaves I (recovery/clearance)
  absorbing <- if (model %in% c("SIS", "SEIS")) states[["S"]] else states[["R"]]

  inf_dist <- .as_dist(infectious_period, "infectious_period")
  inc_dist <- if (has_E) .as_dist(incubation_period, "incubation_period") else NULL
  imm_dist <- if (waning) .as_dist(immunity_period, "immunity_period") else NULL
  beta     <- r0 / mean(inf_dist$sample_n(100000L))

  if (is.null(network)) network <- matrix(0, n_nodes, n_nodes)

  # ── people: state (u8), node id (u16), timer (u16) ─────────────────────────────────
  people <- new.env()
  people$count    <- n
  people$capacity <- n                                  # closed population: no growth
  people$state    <- allocate_scalar("u8",  n)
  people$nodeid   <- allocate_scalar("u16", n)
  people$timer    <- allocate_scalar("u16", n)
  people$nodeid$set(rep(seq_len(n_nodes) - 1L, times = pops))

  I_seed <- if ("I" %in% names(scenario)) as.integer(scenario$I) else integer(n_nodes)
  R_seed <- if ("R" %in% names(scenario)) as.integer(scenario$R) else integer(n_nodes)
  if (any(I_seed + R_seed > pops)) stop("seeded I + R exceeds a node's population")
  offsets <- cumsum(c(0L, pops[-n_nodes]))
  state0  <- rep(states[["S"]], n)
  timer0  <- integer(n)
  for (k in seq_len(n_nodes)) {
    base <- offsets[k]; ni <- I_seed[k]; nr <- R_seed[k]
    if (ni > 0L) {
      state0[base + seq_len(ni)] <- states[["I"]]
      timer0[base + seq_len(ni)] <- .timer0(inf_dist$sample_n(ni))   # infectious clock
    }
    if (nr > 0L) {
      state0[base + ni + seq_len(nr)] <- states[["R"]]
      if (waning) timer0[base + ni + seq_len(nr)] <- .timer0(imm_dist$sample_n(nr))
    }
  }
  people$state$set(state0)
  people$timer$set(timer0)

  # ── nodes: census (only the compartments this model has) + drivers + flows ─────────
  nodes <- new.env()
  nodes$count <- n_nodes
  zeros <- rep(0L, (nticks - 1L) * n_nodes)
  mkcensus <- function(col0) { c <- allocate_vector("i32", nticks, n_nodes); c$set(c(col0, zeros)); c }
  nodes$S <- mkcensus(pops - I_seed - R_seed)
  nodes$I <- mkcensus(I_seed)
  if (has_E) nodes$E <- mkcensus(integer(n_nodes))
  if (has_R) nodes$R <- mkcensus(R_seed)
  nodes$N <- mkcensus(pops)
  nodes$foi       <- allocate_vector("f64", nticks - 1L, n_nodes)
  nodes$incidence <- allocate_vector("i32", nticks - 1L, n_nodes)
  if (has_step_clearance) nodes$recoveries <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$beta        <- values_map(beta,        nticks, n_nodes)
  nodes$seasonality <- values_map(seasonality, nticks, n_nodes)

  # The census compartments to carry forward each tick (and total into N).
  carry <- c(list(nodes$S, nodes$I), if (has_E) list(nodes$E), if (has_R) list(nodes$R))

  # ── per-tick closures ───────────────────────────────────────────────────────────
  progress_step <- function(t) {
    if (model %in% c("SI", "SEI")) {
      r <- step_si(people$state, people$timer, people$nodeid, people$count, n_nodes, inf_dist)
      if (has_E) move_count(nodes$E, nodes$I, r$onset, t)              # E->I
    } else if (waning) {
      r <- step_sirs(people$state, people$timer, people$nodeid, people$count, n_nodes, inf_dist, imm_dist)
      if (has_E) move_count(nodes$E, nodes$I, r$onset, t)              # E->I
      move_count(nodes$I, nodes$R, r$recovered, t)                    # I->R
      move_count(nodes$R, nodes$S, r$waned_r, t)                      # R->S
      nodes$recoveries$set_col(t, r$recovered)
    } else {
      r <- step_sir(people$state, people$timer, people$nodeid, people$count, n_nodes, inf_dist, absorbing)
      if (has_E) move_count(nodes$E, nodes$I, r$onset, t)              # E->I
      to <- if (identical(absorbing, states[["S"]])) nodes$S else nodes$R
      move_count(nodes$I, to, r$cleared, t)                           # I->{S,R}
      nodes$recoveries$set_col(t, r$cleared)
    }
  }
  transmit <- function(t) {
    if (model == "SI") {
      inf <- transmission_si(people$state, people$nodeid, people$count, nodes$foi, t)
      move_count(nodes$S, nodes$I, inf, t)
    } else {
      dest <- if (has_E) states[["E"]] else states[["I"]]
      dist <- if (has_E) inc_dist else inf_dist
      inf <- transmission(people$state, people$timer, people$nodeid, people$count,
                          nodes$foi, t, dest, dist)
      move_count(nodes$S, if (has_E) nodes$E else nodes$I, inf, t)
    }
    nodes$incidence$set_col(t, inf)
  }

  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0L, max = nticks - 1L, style = 3) else NULL
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)
  every <- max(1L, (nticks - 1L) %/% 100L)
  for (tick in seq_len(nticks - 1L)) {
    t <- tick - 1L
    carry_forward_states(carry, t, total = nodes$N)
    if (has_E) progress_step(t)                                       # E-entry: step before FOI
    calc_foi(nodes$I, nodes$N, nodes$beta, nodes$seasonality, network, nodes$foi, t)
    if (!has_E && has_step_clearance) progress_step(t)                # direct-I: step after FOI
    transmit(t)
    if (!is.null(pb) && (tick %% every == 0L || tick == nticks - 1L)) utils::setTxtProgressBar(pb, tick)
  }

  invisible(list(people = people, nodes = nodes))
}
