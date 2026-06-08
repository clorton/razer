# run_model() — a high-level runner for the closed-population compartmental menagerie
# (SI/SEI/SIS/SEIS/SIR/SEIR/SIRS/SEIRS) on the Column kernels. It owns the per-tick
# bookkeeping — the carry-forward, the single `step → calc_foi → transmission` ordering
# that yields R0 = beta*D, and the move_count census deltas — so a user gets a correct
# model from one call instead of ~100 lines of hand-wiring (and there's no way to desync
# the census from the agents).
#
# It bundles `people`, `nodes`, `network`, and `carry` into one `model` environment and
# exposes three optional callbacks — `init(model)` (once, before the loop), and
# `step_enter(model)` / `step_exit(model)` (at the start / end of each tick) — so users can
# add their own properties, compartments, and per-tick logic without forking the runner.
#
# For readers new to R: `switch`/`if` choose per-model behaviour; closures (`function(t)`)
# capture the people/nodes/distribution objects; `inherits(x,"Distribution")` is an
# isinstance check; `grepl("E", model)` tests for a substring; `new.env()` is a mutable,
# reference-semantics record (assigning into it mutates the shared object).

# Coerce a duration argument to a Distribution (a bare number becomes a constant).
.as_dist <- function(x, what) {
  if (inherits(x, "Distribution")) return(x)
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) return(dist_constant(x))
  stop(sprintf("`%s` must be a Distribution or a single finite number", what))
}

# Round duration draws to whole u16 ticks (the kernels' minimum is 1).
.timer0 <- function(x) pmax(1L, pmin(65535L, as.integer(round(x))))

# Validate a count vector is finite, non-negative, and whole-valued; return it as integer.
# Used for the scenario `population` and the E/I/R seed columns so a stray NA, Inf, or
# fractional count fails with a clear message instead of a cryptic base-R error or a silent
# truncation.
.whole_nonneg <- function(x, what) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)))
    stop(sprintf("`%s` must be finite numeric (no NA / Inf)", what))
  if (any(x < 0))         stop(sprintf("`%s` must be non-negative", what))
  if (any(x != floor(x))) stop(sprintf("`%s` must be whole numbers", what))
  as.integer(x)
}

#' Run a closed-population compartmental model.
#'
#' Builds an agent population from a per-node scenario, seeds it, and advances one of the
#' eight SI / SEI / SIS / SEIS / SIR / SEIR / SIRS / SEIRS models for `nticks` ticks using
#' the Column kernels in the correct order (so the realized basic reproduction number is
#' the full `R0 = beta * D`). Each tick runs, uniformly for every model,
#' `carry_forward → step → calc_foi → transmission`, with `calc_foi` placed immediately
#' before `transmission` (no step kernel between them); `run_model` applies every census
#' delta for you.
#'
#' The `people`, `nodes`, `network`, and `carry` objects are bundled into a single `model`
#' environment, which is also what the optional callbacks receive and what is returned.
#'
#' @param scenario A data.frame with one row per node and an integer `population` column.
#'   Optional integer `E` / `I` / `R` columns seed that state per node — but only for
#'   states the model actually has (e.g. an `E` column is ignored by `SIR`). `S` is the
#'   remainder. The seeded `E + I + R` must not exceed a node's population.
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
#' @param capacity Optional agent-array capacity (>= the total initial population, the
#'   default). Reserve extra slots — e.g. `calc_capacity()` / `calc_capacity_cdr()` — so a
#'   callback can grow the population with `births` or `import_infections` (which activate
#'   reserved slots). With the default `capacity`, the population is closed (no growth).
#' @param extra_states Optional character vector of extra agent-state names beyond the
#'   model's own S/E/I/R. A name already in [laser_states()] (e.g. `"M"`) keeps that code; a
#'   NEW name (e.g. `"V"` for vaccinated) is assigned the next free code, becoming a genuine
#'   new agent state. Each gets a census Column that is carried forward, totalled into `N`,
#'   and seeded at tick 0 from the agents' states; the assigned codes are exposed in
#'   `model$states`. The disease kernels branch only on S/E/I/R/M, so an agent in a NEW
#'   state is left untouched (not infected, not transitioned) — its transitions are yours to
#'   drive from callbacks (move `S`→`V` to vaccinate; for waning, set a timer and run
#'   `step_timer_expire(V, S)` in `step_update`). For `"M"` only, run_model additionally
#'   applies the step kernels' built-in M→S waning each tick (recording the `waning_m` flow).
#' @param init Optional `function(model)` called once after the model is built and before
#'   the loop. Use it to add per-agent property [Column]s, add a compartment census Column
#'   (append it to `model$carry` to have it carried forward), set agent states/timers, etc.
#'   The initial per-node census is (re)derived from the agents' states AFTER `init` runs,
#'   so changing agent states in `init` keeps the census consistent.
#' @param step_enter Optional `function(model)` called at the start of every tick (before
#'   the carry-forward). Receives the `model` environment; `model$tick` is the current
#'   0-based interval index.
#' @param step_update Optional `function(model)` called every tick AFTER the step kernel
#'   (the disease-state progression) and BEFORE `calc_foi`. Because `calc_foi` reads the
#'   *settled* start-of-interval census `I[t]` (slice `model$tick`), edit the FOI **drivers**
#'   here — `model$nodes$beta`, `model$nodes$seasonality`, or `model$network` (all read at
#'   slice `t`) — to influence the CURRENT tick's force of infection. Census edits made here
#'   land in the working slice `t+1`, so they affect the NEXT tick's FOI, not this one.
#' @param step_exit Optional `function(model)` called at the end of every tick (after
#'   transmission). Receives the `model` environment (with `model$tick` set).
#' @return Invisibly, the `model` environment, with:
#'   * `$people` — the agent arrays (`state`, `nodeid`, `timer`, plus `count`/`capacity`).
#'   * `$nodes` — per-tick census Columns (`S`, `E`/`I`/`R` as applicable, any `extra_states`
#'     such as `M`), `N`, the `foi`, and the per-interval flows: `incidence` (new infections,
#'     all models), `onset` (E→I, `SE*` models), `recovery` (I to S or R, models with an
#'     I-exit), `waning` (R→S, `SIRS`/`SEIRS`), and `waning_m` (M→S, when `M` is an extra
#'     state) — each an `n_nodes`-wide, time-major Column.
#'   * `$network` — the coupling weights as a 2-D f64 [Column] (`n_nodes x n_nodes`,
#'     column-major), built once from the `network` matrix and read by `calc_foi` each tick.
#'   * `$carry` — the list of census Columns carried forward each tick.
#'   * `$states` — the state-name→code map ([laser_states()] plus any `extra_states`).
#'   * `$tick` — during the loop, the current 0-based interval index (for the callbacks).
#' @examples
#' \dontrun{
#' scenario <- data.frame(population = 1e5L, I = 100L)
#' m <- run_model(scenario, "SEIR", nticks = 200L, r0 = 2.5,
#'                infectious_period = 8, incubation_period = 5, seed = 1)
#' tail(m$nodes$I$values())
#' }
#' @export
run_model <- function(scenario, model, nticks, r0, infectious_period,
                      incubation_period = NULL, immunity_period = NULL,
                      network = NULL, seasonality = 1, seed = NULL, progress = FALSE,
                      capacity = NULL, extra_states = NULL,
                      init = NULL, step_enter = NULL, step_update = NULL, step_exit = NULL) {
  model <- toupper(model)
  valid <- c("SI", "SEI", "SIS", "SEIS", "SIR", "SEIR", "SIRS", "SEIRS")
  if (!model %in% valid)
    stop(sprintf("`model` must be one of %s", paste(valid, collapse = ", ")))
  if (!is.data.frame(scenario) || !"population" %in% names(scenario))
    stop("`scenario` must be a data.frame with an integer `population` column")
  nticks <- as.integer(nticks)
  if (length(nticks) != 1L || is.na(nticks) || nticks < 2L)
    stop("`nticks` must be an integer >= 2")
  if (!is.numeric(r0) || length(r0) != 1L || !is.finite(r0) || r0 < 0)
    stop("`r0` must be a single finite, non-negative number")
  for (cb in c("init", "step_enter", "step_update", "step_exit")) {
    f <- get(cb)
    if (!is.null(f) && !is.function(f)) stop(sprintf("`%s` must be a function or NULL", cb))
  }
  if (!is.null(seed)) set_seed(seed)

  states  <- laser_states()
  pops    <- .whole_nonneg(scenario$population, "scenario$population")
  n       <- sum(pops)
  n_nodes <- nrow(scenario)

  # Agent-array capacity: defaults to the closed-population size `n`. Pass a larger value
  # (e.g. from calc_capacity / calc_capacity_cdr) to reserve slots that a callback can
  # activate with `births` / `import_infections` — i.e. to let the population GROW.
  cap <- if (is.null(capacity)) n else .whole_nonneg(capacity, "capacity")
  if (length(cap) != 1L) stop("`capacity` must be a single number")
  if (cap < n) stop(sprintf("`capacity` (%d) must be >= total initial population (%d)", cap, n))

  has_E   <- grepl("E", model)
  has_R   <- grepl("R", model)
  waning  <- model %in% c("SIRS", "SEIRS")             # I->R sets an immunity timer, R->S
  has_step_clearance <- !model %in% c("SI", "SEI")     # I leaves I (recovery/clearance)
  absorbing <- if (model %in% c("SIS", "SEIS")) states[["S"]] else states[["R"]]

  # Extra compartments beyond the model's own S/E/I/R (e.g. a maternal-immunity "M"). Each
  # is tracked in the node census, carried forward, and totalled into N. For "M" the step
  # kernels already compute the M->S waning leg, so run_model applies it (and records the
  # `waning_m` flow); other extra states are inert until a callback moves agents in/out.
  base_comp <- c("S", "I", if (has_E) "E", if (has_R) "R")
  extra_states <- if (is.null(extra_states)) character(0) else as.character(extra_states)
  if (length(extra_states)) {
    if (any(extra_states %in% c(base_comp, "D")))
      stop("`extra_states` must not repeat the model's own compartments or include D")
    if (anyDuplicated(extra_states)) stop("`extra_states` names must be unique")
    # Give each extra state a u8 code: a known laser_states() name (e.g. "M") keeps its
    # code; a NEW name (e.g. "V" for vaccinated) gets the next free code, becoming a genuine
    # new agent STATE. The disease kernels branch only on S/E/I/R/M, so they leave an agent
    # in a new state UNTOUCHED — not infected, not transitioned, its timer not decremented.
    # Its transitions are yours to drive from callbacks: move S->V to vaccinate; for waning,
    # set a timer on vaccination and run step_timer_expire(V, S) in a step_update callback.
    for (s in extra_states) if (!(s %in% names(states))) {
      free <- setdiff(0:254, as.integer(states))
      if (!length(free)) stop("no free state code available for a new `extra_states` entry")
      states[[s]] <- free[[1L]]
    }
  }
  has_M <- "M" %in% extra_states

  # Warn about inputs the chosen model does not use — usually a typo'd model name or a
  # wrong expectation (e.g. an `E` column or `incubation_period` passed to SIR).
  if (!has_E  && "E" %in% names(scenario))
    warning(sprintf("model %s has no E compartment; the scenario's `E` column is ignored", model))
  if (!has_R  && "R" %in% names(scenario))
    warning(sprintf("model %s has no R compartment; the scenario's `R` column is ignored", model))
  if (!has_E  && !is.null(incubation_period))
    warning(sprintf("model %s has no E compartment; `incubation_period` is ignored", model))
  if (!waning && !is.null(immunity_period))
    warning(sprintf("model %s does not have waning immunity; `immunity_period` is ignored", model))

  inf_dist <- .as_dist(infectious_period, "infectious_period")
  inc_dist <- if (has_E) .as_dist(incubation_period, "incubation_period") else NULL
  imm_dist <- if (waning) .as_dist(immunity_period, "immunity_period") else NULL
  beta     <- r0 / mean(inf_dist$sample_n(100000L))

  if (is.null(network)) network <- matrix(0, n_nodes, n_nodes)
  if (!is.matrix(network) || nrow(network) != n_nodes || ncol(network) != n_nodes)
    stop(sprintf("`network` must be a %d x %d matrix", n_nodes, n_nodes))
  # Copy the coupling matrix into a persistent 2-D f64 Column (column-major, as `as.vector`
  # produces) so calc_foi reads it Rust-side every tick instead of re-marshalling the R
  # matrix. `as.vector` on a matrix yields column-major order, matching calc_foi's indexing.
  network <- {
    nc <- allocate_vector("f64", n_nodes, n_nodes)
    nc$set(as.vector(network))
    nc
  }

  # ── people: state (u8), node id (u16), timer (u16) ─────────────────────────────────
  people <- new.env()
  people$count    <- n
  people$capacity <- cap                                # cap > n leaves reserved growth slots
  people$state    <- allocate_scalar("u8",  cap)
  people$nodeid   <- allocate_scalar("u16", cap)
  people$timer    <- allocate_scalar("u16", cap)
  # Active agents get their node id; reserved slots [n, cap) pad with node 0 (state S, see
  # below) until a birth/import activates them.
  people$nodeid$set(c(rep(seq_len(n_nodes) - 1L, times = pops), integer(cap - n)))

  # ── flexible per-state seeding (item 3) ────────────────────────────────────────────
  # Seed any of the model's E / I / R compartments that the scenario supplies a column
  # for; a state absent from the scenario is left unseeded. S is the remainder. Each
  # seeded state gets its natural timer: E -> incubation, I -> infectious, R -> immunity
  # (only when immunity wanes). States the MODEL lacks are never seeded.
  seed_states <- c(if (has_E) "E", "I", if (has_R) "R")
  seed_counts <- lapply(seed_states, function(st)
    if (st %in% names(scenario)) .whole_nonneg(scenario[[st]], paste0("scenario$", st))
    else integer(n_nodes))
  names(seed_counts) <- seed_states
  total_seeded <- if (length(seed_counts)) Reduce(`+`, seed_counts) else integer(n_nodes)
  if (any(total_seeded > pops))
    stop("seeded E + I + R exceeds a node's population in at least one node")

  offsets <- cumsum(c(0L, pops[-n_nodes]))             # 0-based start of each node's block
  state0  <- rep(states[["S"]], cap)                   # reserved slots default to S / timer 0
  timer0  <- integer(cap)
  for (k in seq_len(n_nodes)) {
    used <- 0L
    for (st in seed_states) {
      cnt <- seed_counts[[st]][k]
      if (cnt > 0L) {
        idx <- offsets[k] + used + seq_len(cnt)
        state0[idx] <- states[[st]]
        if      (st == "I")            timer0[idx] <- .timer0(inf_dist$sample_n(cnt))
        else if (st == "E")            timer0[idx] <- .timer0(inc_dist$sample_n(cnt))
        else if (st == "R" && waning)  timer0[idx] <- .timer0(imm_dist$sample_n(cnt))
        used <- used + cnt
      }
    }
  }
  people$state$set(state0)
  people$timer$set(timer0)

  # ── nodes: census (only this model's compartments) + drivers + flows (item 2) ──────
  nodes <- new.env()
  nodes$count <- n_nodes
  nodes$S <- allocate_vector("i32", nticks, n_nodes)
  nodes$I <- allocate_vector("i32", nticks, n_nodes)
  if (has_E) nodes$E <- allocate_vector("i32", nticks, n_nodes)
  if (has_R) nodes$R <- allocate_vector("i32", nticks, n_nodes)
  for (s in extra_states) nodes[[s]] <- allocate_vector("i32", nticks, n_nodes)  # e.g. M
  nodes$N <- allocate_vector("i32", nticks, n_nodes)
  nodes$foi <- allocate_vector("f64", nticks - 1L, n_nodes)
  nodes$incidence <- allocate_vector("i32", nticks - 1L, n_nodes)              # new infections (all)
  if (has_E)              nodes$onset    <- allocate_vector("i32", nticks - 1L, n_nodes)  # E->I
  if (has_step_clearance) nodes$recovery <- allocate_vector("i32", nticks - 1L, n_nodes)  # I-exit
  if (waning)             nodes$waning   <- allocate_vector("i32", nticks - 1L, n_nodes)  # R->S
  if (has_M)              nodes$waning_m <- allocate_vector("i32", nticks - 1L, n_nodes)  # M->S
  nodes$beta        <- values_map(beta,        nticks, n_nodes)
  nodes$seasonality <- values_map(seasonality, nticks, n_nodes)

  # The census compartments carried forward each tick (and totalled into N). A user can
  # append a compartment Column here in `init` to have it carried too.
  carry <- c(list(nodes$S, nodes$I), if (has_E) list(nodes$E), if (has_R) list(nodes$R),
             lapply(extra_states, function(s) nodes[[s]]))

  # ── bundle into the model environment + run the init callback ──────────────────────
  m <- new.env()
  m$people  <- people
  m$nodes   <- nodes
  m$network <- network
  m$carry   <- carry
  m$states  <- states   # laser_states() plus any extra_states codes — for callbacks
  if (!is.null(init)) init(m)
  # Re-read in case the callback rebound them (it normally mutates in place).
  people <- m$people; nodes <- m$nodes; network <- m$network

  # ── derive the initial (tick 0) census from the agents' states ─────────────────────
  # Done AFTER init so any state changes made there are reflected. `bincount_where` tallies
  # the per-node agent count in each state (and N = the living, state != D == 255).
  nid <- people$nodeid; st <- people$state; cnt <- people$count
  nodes$S$set_col(0L, bincount_where(nid, n_nodes, st, "eq", states[["S"]], cnt))
  nodes$I$set_col(0L, bincount_where(nid, n_nodes, st, "eq", states[["I"]], cnt))
  if (has_E) nodes$E$set_col(0L, bincount_where(nid, n_nodes, st, "eq", states[["E"]], cnt))
  if (has_R) nodes$R$set_col(0L, bincount_where(nid, n_nodes, st, "eq", states[["R"]], cnt))
  for (s in extra_states)
    nodes[[s]]$set_col(0L, bincount_where(nid, n_nodes, st, "eq", states[[s]], cnt))
  nodes$N$set_col(0L, bincount_where(nid, n_nodes, st, "ne", 255, cnt))        # alive (state != D)

  # ── per-tick closures ───────────────────────────────────────────────────────────
  progress_step <- function(t) {
    if (model %in% c("SI", "SEI")) {
      r <- step_si(people$state, people$timer, people$nodeid, people$count, n_nodes, inf_dist)
      if (has_E) { move_count(nodes$E, nodes$I, r$onset, t); nodes$onset$set_col(t, r$onset) }
    } else if (waning) {
      r <- step_sirs(people$state, people$timer, people$nodeid, people$count, n_nodes, inf_dist, imm_dist)
      if (has_E) { move_count(nodes$E, nodes$I, r$onset, t); nodes$onset$set_col(t, r$onset) }
      move_count(nodes$I, nodes$R, r$recovered, t); nodes$recovery$set_col(t, r$recovered)  # I->R
      move_count(nodes$R, nodes$S, r$waned_r, t);   nodes$waning$set_col(t, r$waned_r)       # R->S
    } else {
      r <- step_sir(people$state, people$timer, people$nodeid, people$count, n_nodes, inf_dist, absorbing)
      if (has_E) { move_count(nodes$E, nodes$I, r$onset, t); nodes$onset$set_col(t, r$onset) }
      to <- if (identical(absorbing, states[["S"]])) nodes$S else nodes$R
      move_count(nodes$I, to, r$cleared, t); nodes$recovery$set_col(t, r$cleared)            # I->{S,R}
    }
    # All step kernels lead with the maternal M->S leg; apply it when an M compartment is
    # registered (run_model otherwise discards `waned`, since the menagerie has no M).
    if (has_M) { move_count(nodes$M, nodes$S, r$waned, t); nodes$waning_m$set_col(t, r$waned) }
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
    m$tick <- t                                                       # 0-based interval index for callbacks
    if (!is.null(step_enter)) step_enter(m)
    carry_forward_states(m$carry, t, total = nodes$N)
    progress_step(t)                                                  # step always before FOI
    if (!is.null(step_update)) step_update(m)                         # between step and FOI
    calc_foi(nodes$I, nodes$N, nodes$beta, nodes$seasonality, network, nodes$foi, t)
    transmit(t)                                                       # FOI immediately before transmit
    if (!is.null(step_exit)) step_exit(m)
    if (!is.null(pb) && (tick %% every == 0L || tick == nticks - 1L)) utils::setTxtProgressBar(pb, tick)
  }

  invisible(m)
}
