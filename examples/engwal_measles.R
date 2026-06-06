#!/usr/bin/env Rscript

# England & Wales measles model (work in progress). Measles needs a richer compartment
# structure than the SIR/SEIR examples:
#
#   M  maternal immunity  newborns temporarily protected by maternal antibodies; after a
#                         waning period (~9 months) they become Susceptible.
#   S  susceptible
#   E  exposed/incubating infected but not yet infectious; after a latent period -> I.
#   I  infectious
#   R  recovered          lifelong immunity.
#   D  deceased           absorbing; removed from the living compartments.
#
# Disease progression is M -> S -> E -> I -> R, plus vital dynamics (a crude birth rate
# feeds newborns into M) and natural (non-disease) mortality. THIS revision adds the
# natural-mortality machinery: every agent is given a realistic age (from an age
# distribution curve) and a date of death `dod` (from a Kaplan-Meier estimator built on
# the SAME curve); a per-tick mortality() kernel retires agents when their `dod` arrives.
#
# Run from anywhere:  Rscript examples/engwal_measles.R

library(razer)

# ── run_measles_model: assemble parameters and run the dynamics ──────────────────────
# Like run_sir_model() (examples/simple_sir.R). Parameters:
#   scenario             data.frame, one row per patch (name, population), with optional
#                        integer seed columns `I` / `R`.
#   network              N x N spatial coupling matrix (1x1 zeros for a single patch).
#   nticks               number of recorded daily states (dynamics run nticks-1 times).
#   inf_duration         Distribution for the infectious period (I -> R timer).
#   incubation_duration  Distribution for the latent/exposed period (E -> I timer).
#   maternal_duration    Distribution for maternal-immunity waning (M -> S timer); needs
#                        a uint16 timer (~270 days > a uint8's 255 range).
#   r0                   basic reproduction number; beta = r0 / mean(inf_duration).
#   cbr                  crude birth rate (births per 1,000 per year); sizes capacity and
#                        (later) drives births into M.
#   age_dist             AliasedDistribution over single-year age bins (initial ages).
#   km                   KaplanMeierEstimator for sampling each agent's age at death.
#   progress             draw a text progress bar over the per-tick loop.
#
# Per-tick structure (measles = SEIR + M, so step_sir with absorbing = R). Each kernel
# returns per-node counts that the model applies to the census with move_count():
#   carry_forward_states(list(M, S, E, I, R), t0, total = N)
#   step_sir(absorbing = R) # M->S waning, E->I incubation, I->R recovery (one u16 pass)
#   calc_foi(...)           # force of infection from the settled census I[t0] / N[t0]
#   transmission(... E)     # S->E exposures, drawing the incubation timer
#   births(...)             # CBR newborns into M (maternal timer + a KM date of death)
#   mortality(...)          # retire agents whose dod has arrived (deaths by compartment)
#
# Returns (invisibly) a list with the `people` and `nodes` environments.
run_measles_model <- function(scenario, network, nticks,
                              inf_duration, incubation_duration, maternal_duration,
                              r0, cbr, age_dist, km, progress = FALSE) {
  if (!inherits(inf_duration, "Distribution") ||
      !inherits(incubation_duration, "Distribution") ||
      !inherits(maternal_duration, "Distribution"))
    stop("`inf_duration`, `incubation_duration`, `maternal_duration` must be Distributions")
  if (!inherits(age_dist, "AliasedDistribution"))
    stop("`age_dist` must be an AliasedDistribution (from aliased_distribution())")
  if (!inherits(km, "KaplanMeierEstimator"))
    stop("`km` must be a KaplanMeierEstimator (from kaplan_meier_estimator())")

  states  <- laser_states()              # c(S=0, E=1, I=2, R=3, M=4, D=-1)
  pops    <- scenario$population
  n       <- sum(pops)                   # initial live-agent count (sum over patches)
  n_nodes <- nrow(scenario)

  # ── capacity for the full run ─────────────────────────────────────────────────────
  # Births grow the population, so allocate room for every agent that will ever exist.
  # calc_capacity() projects births-driven growth from the (constant) CBR over all
  # nticks-1 daily steps; sum the per-node estimates for the total slots to allocate.
  birthrates <- matrix(cbr, nrow = nticks - 1L, ncol = n_nodes)
  capacity   <- as.integer(sum(calc_capacity(birthrates, pops, safety_factor = 1)))

  # Transmission coefficient: R0 = beta * mean infectious period, so beta = R0 / mean.
  # (calc_foi reads the settled start-of-interval census I[t0], so each infectious agent
  # contributes on exactly the D census columns it occupies — the full infectious period
  # D counts and R0 = beta * D — see CLAUDE.md.)
  beta <- r0 / mean(inf_duration$sample_n(100000L))

  # ── people: per-agent property arrays, allocated to capacity ───────────────────────
  # state (u8), nodeid (u16), a uint16 `timer` (maternal waning needs > 255 ticks), and
  # the demographic clocks: `dob` (date of birth, a SIGNED i32 — negative for agents born
  # before the simulation starts) and `dod` (date of death, an absolute u32 tick).
  people <- new.env()
  people$count    <- n
  people$capacity <- capacity
  people$state    <- allocate_scalar("u8",  capacity)
  people$nodeid   <- allocate_scalar("u16", capacity)
  people$timer    <- allocate_scalar("u16", capacity)
  people$dob      <- allocate_scalar("i32", capacity)
  people$dod      <- allocate_scalar("u32", capacity)

  # nodeid is 0-based; repeat each node id by its population, padding reserved slots with 0.
  ids <- seq_len(n_nodes) - 1L
  people$nodeid$set(c(rep(ids, times = pops), rep(0L, capacity - n)))

  # ── ages, dates of birth, and dates of death ───────────────────────────────────────
  # Sample each initial agent's age (in days): a 0-based age-year bin from `age_dist`,
  # then a uniform day within that year. Date of birth = -age (born `age` days ago). Age
  # AT death comes from the Kaplan-Meier estimator (conditioned on the current age);
  # `dod` is that minus the current age = the absolute tick at which the agent will die.
  age_days     <- age_dist$sample_n(n) * 365L + as.integer(floor(stats::runif(n) * 365))
  age_at_death <- km$predict_age_at_death(age_days, -1L)   # -1L: use the full life table
  dob_vals     <- -age_days
  dod_vals     <- age_at_death - age_days                  # >= 0 (death not in the past)
  people$dob$set(c(dob_vals, rep(0L,  capacity - n)))
  people$dod$set(c(dod_vals, rep(0L,  capacity - n)))

  # ── seed the initial disease states: I seed, R immune fraction, the rest S ─────────
  # Agents are laid out node-by-node; node k owns [offset[k], offset[k] + pops[k]). M and
  # E start empty (newborns will fill M as births are added later).
  I_seed  <- if ("I" %in% names(scenario)) as.integer(scenario$I) else integer(n_nodes)
  R_seed  <- if ("R" %in% names(scenario)) as.integer(scenario$R) else integer(n_nodes)
  offsets <- cumsum(c(0L, pops[-n_nodes]))
  state0  <- rep(states[["S"]], capacity)
  for (k in seq_len(n_nodes)) {
    base <- offsets[k]; ni <- I_seed[k]; nr <- R_seed[k]
    if (ni > 0L) state0[base + seq_len(ni)]      <- states[["I"]]
    if (nr > 0L) state0[base + ni + seq_len(nr)] <- states[["R"]]
  }
  people$state$set(state0)

  # Seed the infectious timer (u16) for the initial I agents from `inf_duration`, so
  # step_sir counts them down to recovery rather than recovering them immediately
  # (a timer of 0 would expire on the first step). Everyone else starts at 0 (M/E start
  # empty; S and R are untimed).
  timer0  <- rep(0L, capacity)
  total_I <- sum(I_seed)
  if (total_I > 0L) {
    draws <- pmax(1L, pmin(65535L, as.integer(round(inf_duration$sample_n(total_I)))))
    pos <- 0L
    for (k in seq_len(n_nodes)) {
      ni <- I_seed[k]
      if (ni > 0L) { timer0[offsets[k] + seq_len(ni)] <- draws[pos + seq_len(ni)]; pos <- pos + ni }
    }
  }
  people$timer$set(timer0)

  # ── nodes: census (M/S/E/I/R + N) and the deaths flow ──────────────────────────────
  nodes <- new.env()
  nodes$count  <- n_nodes
  nodes$M      <- allocate_vector("i32", nticks,      n_nodes)
  nodes$S      <- allocate_vector("i32", nticks,      n_nodes)
  nodes$E      <- allocate_vector("i32", nticks,      n_nodes)
  nodes$I      <- allocate_vector("i32", nticks,      n_nodes)
  nodes$R      <- allocate_vector("i32", nticks,      n_nodes)
  nodes$N      <- allocate_vector("i32", nticks,      n_nodes)
  # Per-interval FLOW reports and the FOI working buffer.
  nodes$foi       <- allocate_vector("f64", nticks - 1L, n_nodes)
  nodes$incidence <- allocate_vector("i32", nticks - 1L, n_nodes)   # S->E exposures
  nodes$births    <- allocate_vector("i32", nticks - 1L, n_nodes)
  nodes$deaths    <- allocate_vector("i32", nticks - 1L, n_nodes)
  # Transmission driver grids (n_ticks x n_nodes f64), built by values_map: the beta
  # coefficient, a (here flat) seasonal modifier, and the daily birth rate (CBR per
  # 1,000 per year -> per person per day) feeding births.
  nodes$beta        <- values_map(beta,             nticks, n_nodes)
  nodes$seasonality <- values_map(1,                nticks, n_nodes)
  nodes$birth_rate  <- values_map(cbr / 1000 / 365, nticks, n_nodes)

  # Seed the census at tick 0 (first n_nodes entries); remaining columns start zero.
  zeros <- rep(0L, (nticks - 1L) * n_nodes)
  nodes$M$set(c(integer(n_nodes),         zeros))   # M_0 = 0
  nodes$S$set(c(pops - I_seed - R_seed,   zeros))
  nodes$E$set(c(integer(n_nodes),         zeros))   # E_0 = 0
  nodes$I$set(c(I_seed,                   zeros))
  nodes$R$set(c(R_seed,                   zeros))
  nodes$N$set(c(pops,                     zeros))

  # ── per-tick loop ───────────────────────────────────────────────────────────────
  run <- function() {
    pb <- if (progress) utils::txtProgressBar(min = 0L, max = nticks - 1L, style = 3) else NULL
    on.exit(if (!is.null(pb)) close(pb), add = TRUE)
    update_every <- max(1L, (nticks - 1L) %/% 100L)
    for (tick in seq_len(nticks - 1L)) {
      t0 <- tick - 1L
      # Carry M/S/E/I/R forward and total them into N (the FOI denominator).
      carry_forward_states(list(nodes$M, nodes$S, nodes$E, nodes$I, nodes$R), t0, total = nodes$N)
      # Timed transitions in one pass (measles = SEIR + M, so step_sir absorbing = R):
      # M->S (waning), E->I (incubation; draws an infectious timer), I->R (recovery).
      # The kernel returns per-node counts; apply each to the census.
      prog <- step_sir(people$state, people$timer, people$nodeid, people$count,
                       nodes$count, inf_duration, states[["R"]])
      move_count(nodes$M, nodes$S, prog$waned,   t0)   # M -> S
      move_count(nodes$E, nodes$I, prog$onset,   t0)   # E -> I
      move_count(nodes$I, nodes$R, prog$cleared, t0)   # I -> R
      # Force of infection from the SETTLED start-of-interval census I[t0]/N[t0], placed
      # IMMEDIATELY before transmission, then S->E exposures (drawing the incubation
      # timer); incidence records new exposures. Reading the settled census gives each
      # infectious agent its full D census columns, so R0 = beta * D.
      calc_foi(nodes$I, nodes$N, nodes$beta, nodes$seasonality, network, nodes$foi, t0)
      inf <- transmission(people$state, people$timer, people$nodeid, people$count,
                          nodes$foi, t0, states[["E"]], incubation_duration)
      move_count(nodes$S, nodes$E, inf, t0)
      nodes$incidence$set_col(t0, inf)
      # Births: CBR newborns into M (maternal timer + a Kaplan-Meier date of death);
      # `births` returns the grown count and the per-node birth count.
      b <- births(people$state, people$timer, people$nodeid, people$dob, people$dod,
                  people$count, nodes$count, nodes$birth_rate, maternal_duration, km, t0)
      people$count <- b$count
      move_count(NULL, nodes$M, b$born, t0)            # newborns enter M
      nodes$births$set_col(t0, b$born)
      # Natural mortality: returns deaths per node by source compartment; decrement each.
      d <- mortality(people$state, people$dod, people$nodeid, people$count, nodes$count, t0)
      move_count(nodes$M, NULL, d$m, t0); move_count(nodes$S, NULL, d$s, t0)
      move_count(nodes$E, NULL, d$e, t0); move_count(nodes$I, NULL, d$i, t0)
      move_count(nodes$R, NULL, d$r, t0)
      nodes$deaths$set_col(t0, d$m + d$s + d$e + d$i + d$r)
      if (!is.null(pb) && (tick %% update_every == 0L || tick == nticks - 1L))
        utils::setTxtProgressBar(pb, tick)
    }
  }
  run()

  # Report. Living population per tick = M+S+E+I+R; final is row `nticks`.
  living <- function(t) sum(nodes$M$values()[t, ], nodes$S$values()[t, ], nodes$E$values()[t, ],
                            nodes$I$values()[t, ], nodes$R$values()[t, ])
  cat(sprintf(paste0("run_measles_model: %d node(s), %s initial agents (capacity %s), %d ticks (%.0f yr); ",
                     "living %s -> %s; cases %s, births %s, natural deaths %s.\n"),
              n_nodes, format(n, big.mark = ","), format(capacity, big.mark = ","),
              nticks, nticks / 365,
              format(living(1L), big.mark = ","), format(living(nticks), big.mark = ","),
              format(sum(nodes$incidence$values()), big.mark = ","),
              format(sum(nodes$births$values()), big.mark = ","),
              format(sum(nodes$deaths$values()), big.mark = ",")))
  invisible(list(people = people, nodes = nodes))
}

# ── setup ───────────────────────────────────────────────────────────────────────────
# Resolve this script's directory so plots land in examples/output/ regardless of the
# working directory.
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# One node of a million agents.
scenario <- data.frame(name = "engwal", population = 1000000L)
n_nodes  <- nrow(scenario)
network  <- matrix(0, nrow = n_nodes, ncol = n_nodes)   # single patch: no coupling

r0                  <- 14
inf_duration        <- dist_normal(7, 1.5)     # infectious period (mean 7, variance 1.5)
incubation_duration <- dist_normal(7, 1.5)     # latent/exposed period
maternal_duration   <- dist_normal(270, 20)    # maternal-immunity waning (~9 months)
cbr                 <- 30                       # crude birth rate (per 1,000 per year)
nticks              <- 20L * 365L               # twenty years of daily steps

# ── age distribution curve + life table (the SAME curve drives both) ─────────────────
# Build a survivorship curve l(a) from a Gompertz-Makeham per-year mortality hazard (a
# small baseline + exponential rise + elevated infant rate). The number alive at each age
# is the age-distribution curve we sample initial ages from; its complement,
# cumulative_deaths(a) = l(0) - l(a+1), is the life table the Kaplan-Meier estimator uses
# to draw dates of death. Deriving both from one curve keeps the population's age
# structure and its mortality mutually consistent.
max_age   <- 100L
age_years <- 0:max_age
hazard    <- 0.0004 + 1e-5 * exp(0.115 * age_years)   # Gompertz-Makeham-ish
hazard[1] <- 0.02                                     # infant (age 0) mortality
hazard    <- pmin(hazard, 1)

surv <- numeric(max_age + 2L)                         # l(0 .. max_age+1)
surv[1] <- 1
for (a in seq_len(max_age + 1L)) surv[a + 1L] <- surv[a] * (1 - hazard[a])
surv[max_age + 2L] <- 0                               # force everyone dead by the top age

cohort            <- 1e6                              # scale fractions to integer counts
age_curve         <- surv[1:(max_age + 1L)]           # l(0..max_age): the age distribution
age_counts        <- round(age_curve * cohort)        # weights for the alias sampler
cumulative_deaths <- round((surv[1] - surv[2:(max_age + 2L)]) * cohort)  # life table

age_dist <- aliased_distribution(age_counts)
km       <- kaplan_meier_estimator(cumulative_deaths)

# Plot the age distribution curve and its life table for review.
grDevices::png(file.path(out_dir, "engwal_measles_age_curve.png"),
               width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))
plot(age_years, age_counts / sum(age_counts), type = "h", lwd = 3, col = "steelblue",
     xlab = "age (years)", ylab = "fraction of population",
     main = "Age distribution curve (initial-age sampling weights)")
plot(age_years, cumulative_deaths / cohort, type = "l", lwd = 2, col = "firebrick",
     xlab = "age (years)", ylab = "cumulative fraction dead",
     main = "Implied life table (Kaplan-Meier source: cumulative deaths by age)")
graphics::par(op)
grDevices::dev.off()
cat(sprintf("wrote age-distribution / life-table plot to %s\n",
            file.path(out_dir, "engwal_measles_age_curve.png")))

# Initial conditions: immune fraction 1 - 1/R0 (herd-immunity threshold), small I seed.
scenario$R <- round((1 - 1 / r0) * scenario$population)
scenario$I <- 100L

# ── run ───────────────────────────────────────────────────────────────────────────
timing <- system.time(
  result <- run_measles_model(scenario, network, nticks,
                              inf_duration, incubation_duration, maternal_duration,
                              r0, cbr, age_dist, km, progress = TRUE))
cat(sprintf("run_measles_model completed in %.3f s\n", timing[["elapsed"]]))

# ── plot the measles dynamics ───────────────────────────────────────────────────────
# `$values()` returns each census buffer as an n_ticks x n_nodes matrix; rowSums totals
# over patches. Plot the large immune-status reservoirs (S, R, M) and the epidemic curve
# (daily incidence = new S->E exposures) on a shared time axis in years.
years           <- (seq_len(nticks) - 1L) / 365
S <- rowSums(result$nodes$S$values()); R <- rowSums(result$nodes$R$values())
M <- rowSums(result$nodes$M$values()); I <- rowSums(result$nodes$I$values())
incidence_daily <- rowSums(result$nodes$incidence$values())

grDevices::png(file.path(out_dir, "engwal_measles_dynamics.png"),
               width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))

# Top: the immune-status reservoirs. `matplot` draws each column as its own line.
matplot(years, cbind(S, R, M), type = "l", lty = 1, lwd = 2,
        col = c("steelblue", "darkgreen", "purple"),
        xlab = "time (years)", ylab = "agents",
        main = "Measles reservoirs: susceptible, recovered, maternally-immune")
graphics::abline(h = sum(scenario$population) / r0, col = "steelblue", lty = 3)  # S* = N/R0
legend("right", legend = c("S", "R", "M", "N / R0"), bty = "n",
       col = c("steelblue", "darkgreen", "purple", "steelblue"),
       lwd = c(2, 2, 2, 1), lty = c(1, 1, 1, 3))

# Bottom: the epidemic curve (daily new infections), with year gridlines.
plot(years[-1], incidence_daily, type = "l", col = "firebrick",
     xlab = "time (years)", ylab = "new infections / day",
     main = "Daily incidence (S -> E exposures)")
graphics::abline(v = 0:20, col = "grey90")
graphics::par(op)
grDevices::dev.off()
cat(sprintf("wrote measles dynamics plot to %s\n",
            file.path(out_dir, "engwal_measles_dynamics.png")))
