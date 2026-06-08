#!/usr/bin/env Rscript

# England & Wales measles model. Measles needs a richer compartment structure than the
# SIR/SEIR examples:
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
# feeds newborns into M) and natural (non-disease) mortality (every agent has a realistic
# age and a Kaplan-Meier date of death; a mortality() kernel retires them when it arrives).
#
# Built on the high-level runner run_model(), exercising the two features that let it go
# beyond the closed-population menagerie:
#   * `extra_states = "M"` — registers the maternal compartment. run_model tracks M in the
#     census, carries it forward, totals it into N, and applies the step kernels' built-in
#     M -> S waning each tick (so a newborn is protected until its maternal timer expires).
#   * `capacity` — sized with calc_capacity() so the births kernel has reserved slots to
#     grow the population into over the 20-year run.
# The disease transitions (M->S, E->I, I->R, S->E) are the runner's; births-into-M and
# natural mortality are added through a step_exit callback.
#
# Run from anywhere:  Rscript examples/engwal_measles.R

library(razer)

# ── run_measles_model: assemble parameters and run via run_model() ────────────────────
# `scenario` data.frame (name, population, optional I / R seed columns); `network` the
# N x N coupling matrix (1x1 zeros for a single patch); `nticks` recorded daily states;
# `inf_duration` / `incubation_duration` / `maternal_duration` the period Distributions;
# `r0` the basic reproduction number; `cbr` the crude birth rate (per 1,000/year);
# `age_dist` an AliasedDistribution over age bins; `km` a KaplanMeierEstimator. Returns
# (invisibly) the `model` environment.
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

  pops    <- scenario$population
  n       <- sum(pops)
  n_nodes <- nrow(scenario)

  # Births grow the population, so size the agent arrays for everyone ever born over the
  # run: calc_capacity projects births-driven growth from the (constant) CBR. run_model
  # reserves these slots; the births kernel activates them.
  birthrates <- matrix(cbr, nrow = nticks - 1L, ncol = n_nodes)
  capacity   <- as.integer(sum(calc_capacity(birthrates, pops, safety_factor = 1)))

  # SEIR + a maternal M compartment. run_model owns the M->S / E->I / I->R disease step and
  # the S->E transmission (drawing the incubation timer); M->S is applied because "M" is a
  # registered extra state. Births-into-M and natural mortality are the step_exit callback.
  m <- run_model(
    scenario = scenario, model = "SEIR", nticks = nticks, r0 = r0,
    infectious_period = inf_duration, incubation_period = incubation_duration,
    network = network, capacity = capacity, extra_states = "M", seed = 1L, progress = progress,
    init = function(model) {
      cap <- model$people$capacity
      model$people$dob <- allocate_scalar("i32", cap)   # date of birth (negative = before t0)
      model$people$dod <- allocate_scalar("u32", cap)   # absolute tick of death
      model$nodes$birth_rate <- values_map(cbr / 1000 / 365, nticks, n_nodes)
      model$nodes$births     <- allocate_vector("i32", nticks - 1L, n_nodes)
      model$nodes$deaths     <- allocate_vector("i32", nticks - 1L, n_nodes)
      # Realistic initial ages from the age curve, and a KM age-at-death per agent
      # (conditioned on the current age, so nobody dies in the past). dob = -age; dod is the
      # absolute tick of death. Reserved slots [n, cap) stay 0 until a birth fills them.
      age_days     <- age_dist$sample_n(n) * 365L + as.integer(floor(stats::runif(n) * 365))
      age_at_death <- km$predict_age_at_death(age_days, -1L)
      model$people$dob$set(c(-age_days, integer(cap - n)))
      model$people$dod$set(c(as.integer(age_at_death - age_days), integer(cap - n)))
    },
    step_exit = function(model) {
      t <- model$tick
      # Births: CBR newborns into M (a maternal timer + a Kaplan-Meier date of death).
      b <- births(model$people$state, model$people$timer, model$people$nodeid,
                  model$people$dob, model$people$dod, model$people$count, n_nodes,
                  model$nodes$birth_rate, maternal_duration, km, t)
      model$people$count <- b$count
      move_count(NULL, model$nodes$M, b$born, t)
      model$nodes$births$set_col(t, b$born)
      # Natural mortality: retire agents whose dod has arrived; decrement their compartments.
      d <- mortality(model$people$state, model$people$dod, model$people$nodeid,
                     model$people$count, n_nodes, t)
      move_count(model$nodes$M, NULL, d$m, t); move_count(model$nodes$S, NULL, d$s, t)
      move_count(model$nodes$E, NULL, d$e, t); move_count(model$nodes$I, NULL, d$i, t)
      move_count(model$nodes$R, NULL, d$r, t)
      model$nodes$deaths$set_col(t, d$m + d$s + d$e + d$i + d$r)
    })

  # Report. Living population per tick = M+S+E+I+R; final is row `nticks`.
  living <- function(t) sum(m$nodes$M$values()[t, ], m$nodes$S$values()[t, ], m$nodes$E$values()[t, ],
                            m$nodes$I$values()[t, ], m$nodes$R$values()[t, ])
  cat(sprintf(paste0("run_measles_model: %d node(s), %s initial agents (capacity %s), %d ticks (%.0f yr); ",
                     "living %s -> %s; cases %s, births %s, natural deaths %s.\n"),
              n_nodes, format(n, big.mark = ","), format(capacity, big.mark = ","),
              nticks, nticks / 365,
              format(living(1L), big.mark = ","), format(living(nticks), big.mark = ","),
              format(sum(m$nodes$incidence$values()), big.mark = ","),
              format(sum(m$nodes$births$values()), big.mark = ","),
              format(sum(m$nodes$deaths$values()), big.mark = ",")))
  invisible(m)
}

# ── setup ───────────────────────────────────────────────────────────────────────────
# Resolve this script's directory so plots land in examples/output/ regardless of the
# working directory.
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."
# Device-aware: write a PNG when run non-interactively (Rscript); draw to the active
# device (e.g. RStudio's Plots pane) when sourced interactively.
to_png    <- !interactive()
open_png  <- function(path, ...) if (to_png) grDevices::png(path, ...)
close_png <- function() if (to_png) grDevices::dev.off()
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
open_png(file.path(out_dir, "engwal_measles_age_curve.png"),
               width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))
plot(age_years, age_counts / sum(age_counts), type = "h", lwd = 3, col = "steelblue",
     xlab = "age (years)", ylab = "fraction of population",
     main = "Age distribution curve (initial-age sampling weights)")
plot(age_years, cumulative_deaths / cohort, type = "l", lwd = 2, col = "firebrick",
     xlab = "age (years)", ylab = "cumulative fraction dead",
     main = "Implied life table (Kaplan-Meier source: cumulative deaths by age)")
graphics::par(op)
close_png()
if (to_png) cat(sprintf("wrote age-distribution / life-table plot to %s\n",
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

open_png(file.path(out_dir, "engwal_measles_dynamics.png"),
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
close_png()
if (to_png) cat(sprintf("wrote measles dynamics plot to %s\n",
            file.path(out_dir, "engwal_measles_dynamics.png")))
