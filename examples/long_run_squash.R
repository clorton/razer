#!/usr/bin/env Rscript

# A century-long demographic run that stays within a bounded agent array by reclaiming the
# slots of the dead with squash(). This is the use case calc_capacity_cdr() exists for.
#
# Setup: one node, 1,000,000 initial agents, 100 years (daily steps), crude BIRTH rate
# 30 / 1,000 / year and crude DEATH rate 15 / 1,000 / year (so the population grows ~15 /
# 1,000 / year net). With no reclaim you would allocate a slot for every agent EVER born —
# tens of millions over a century (what calc_capacity() estimates). Instead we size the
# array with calc_capacity_cdr() (the PEAK-LIVING bound) and squash() once a year to free
# the dead slots, so the active count stays near the living count, well under capacity.
#
# Built on run_model(): there is no disease here, so we run the trivial model "SI" with
# beta = 0 (no transmission). The vital dynamics are added through run_model's callbacks:
#   * `capacity` reserves the (calc_capacity_cdr) headroom births need to grow into;
#   * `extra_states = "M"` registers the maternal-immunity state newborns enter
#     (run_model applies the step kernel's M->S waning each tick);
#   * `init` gives each agent a date of birth / Kaplan-Meier date of death;
#   * `step_exit` runs births (into M) and mortality each tick, and squash() once a year.
#
# Run from anywhere:  Rscript examples/long_run_squash.R   (takes a few minutes)
# Output PNGs are written next to this script in examples/output/.

library(razer)

set.seed(20260606L)

# ── parameters ────────────────────────────────────────────────────────────────────
N0        <- 1000000L                     # initial living population
years     <- 100L
nticks    <- years * 365L                 # daily steps
cbr       <- 30                           # crude birth rate (per 1,000 per year)
cdr       <- 15                           # crude death rate (per 1,000 per year)
maternal_duration <- dist_normal(270, 20) # M -> S waning (~9 months)
squash_every <- 365L                      # reclaim dead slots once a year

# ── capacity: peak-living bound (with squash) vs cumulative-births bound (without) ──
br_grid   <- matrix(cbr, nrow = nticks - 1L, ncol = 1L)
dr_grid   <- matrix(cdr, nrow = nticks - 1L, ncol = 1L)
capacity  <- as.integer(calc_capacity_cdr(br_grid, dr_grid, N0, safety_factor = 1))
cap_naive <- as.integer(min(calc_capacity(br_grid, N0, safety_factor = 1), .Machine$integer.max))
cat(sprintf("capacity (squash-aware calc_capacity_cdr) = %s slots\n", format(capacity, big.mark = ",")))
cat(sprintf("capacity (no-reclaim calc_capacity)       = %s slots  (%.1fx more)\n",
            format(cap_naive, big.mark = ","), cap_naive / capacity))

# ── exponential life table for crude (age-independent) mortality at the CDR ─────────
# Constant hazard h = -log(1 - CDR/1000) per year gives an exponential lifespan, so the
# per-capita death rate is ~CDR at every age. The KM estimator samples ages at death from
# the cumulative-deaths curve (used for both initial agents and newborns).
max_age   <- 200L
h_year    <- -log(1 - cdr / 1000)
surv      <- exp(-h_year * (0:max_age))           # l(a) = exp(-h a)
cohort    <- 1e7
cumulative_deaths <- round((1 - surv) * cohort)
km        <- kaplan_meier_estimator(cumulative_deaths)

# ── per-tick series recorded by the step_exit callback (env so the closure can write) ─
# `track$live[i]` / `track$active[i]` hold the living population and the active agent count
# (people$count) at census slice i-1; `births` / `deaths` are per-interval flows.
track <- new.env()
track$live   <- numeric(nticks); track$live[1]   <- N0
track$active <- numeric(nticks); track$active[1] <- N0
track$births <- numeric(nticks - 1L)
track$deaths <- numeric(nticks - 1L)

# ── run the 100-year demographic model through run_model() ──────────────────────────
timing <- system.time(
  model <- run_model(
    scenario = data.frame(population = N0), model = "SI", nticks = nticks, beta = 0,
    infectious_period = dist_constant(7),     # unused: beta = 0 means no transmission
    capacity = capacity, extra_states = "M", seed = 1L, progress = TRUE,
    init = function(m) {
      cap <- m$people$capacity
      m$people$dob <- allocate_scalar("i32", cap)        # date of birth (negative = before t0)
      m$people$dod <- allocate_scalar("u32", cap)        # absolute tick of death
      m$nodes$birth_rate <- values_map(cbr / 1000 / 365, nticks, m$nodes$count)
      # Initial ages from the exponential stationary structure (capped at 90 years, beyond
      # the life table's reach); dod = now + a memoryless remaining lifespan from the table.
      age_days     <- pmin(as.integer(floor(stats::rexp(N0, rate = h_year / 365))), 90L * 365L)
      age_at_death <- km$predict_age_at_death(age_days, -1L)
      m$people$dob$set(c(-age_days, integer(cap - N0)))
      m$people$dod$set(c(as.integer(age_at_death - age_days), integer(cap - N0)))
    },
    step_exit = function(m) {
      t <- m$tick; bi <- t + 1L; idx <- t + 2L
      # Births: CBR newborns into M (maternal timer + a KM date of death). run_model wanes
      # M -> S each tick (the "M" extra state), so we only ADD newborns here.
      b <- births(m$people$state, m$people$timer, m$people$nodeid, m$people$dob, m$people$dod,
                  m$people$count, m$nodes$count, m$nodes$birth_rate, maternal_duration, km, t)
      m$people$count <- b$count
      move_count(NULL, m$nodes$M, b$born, t)
      track$births[bi] <- sum(b$born)
      # Crude mortality: retire agents whose dod has arrived; decrement their states
      # (only S, M, and the unused I are populated here).
      d <- mortality(m$people$state, m$people$dod, m$people$nodeid, m$people$count, m$nodes$count, t)
      move_count(m$nodes$M, NULL, d$m, t)
      move_count(m$nodes$S, NULL, d$s, t)
      move_count(m$nodes$I, NULL, d$i, t)
      track$deaths[bi] <- sum(d$m + d$s + d$i)
      # Once a year, reclaim the slots of the dead so the active count tracks the LIVING
      # count instead of growing toward the cumulative-births total.
      if (bi %% squash_every == 0L) m$people$count <- squash(m$people)
      track$live[idx]   <- sum(m$nodes$M$values()[idx, ], m$nodes$S$values()[idx, ])
      track$active[idx] <- m$people$count
    }))
cat(sprintf("\nrun took %.1f s\n", timing[["elapsed"]]))

# ── report ────────────────────────────────────────────────────────────────────────
live_pop <- track$live; active <- track$active; births_d <- track$births; deaths_d <- track$deaths
person_years <- sum(live_pop) / 365
cat(sprintf("living population: %s -> %s over %d years (net growth %.2fx)\n",
            format(N0, big.mark = ","), format(round(live_pop[nticks]), big.mark = ","),
            years, live_pop[nticks] / N0))
cat(sprintf("total births = %s, total deaths = %s\n",
            format(sum(births_d), big.mark = ","), format(sum(deaths_d), big.mark = ",")))
cat(sprintf("realized CBR = %.1f, CDR = %.1f per 1,000/yr (requested %d / %d)\n",
            1000 * sum(births_d) / person_years, 1000 * sum(deaths_d) / person_years, cbr, cdr))
cat(sprintf("peak active slots = %s of %s capacity (%.0f%% utilization); peak living = %s\n",
            format(round(max(active)), big.mark = ","), format(capacity, big.mark = ","),
            100 * max(active) / capacity, format(round(max(live_pop)), big.mark = ",")))

# ── output location ─────────────────────────────────────────────────────────────────
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

yrs <- (seq_len(nticks) - 1L) / 365

# ── plot 1: living population vs the two capacity bounds ─────────────────────────────
open_png(file.path(out_dir, "long_run_squash_population.png"), width = 1050, height = 680, res = 110)
plot(yrs, live_pop / 1e6, type = "l", lwd = 2.5, col = "#238b45",
     ylim = c(0, cap_naive / 1e6 * 1.05), xlab = "year", ylab = "agents (millions)",
     main = "100-year run: living population vs. agent-array capacity bounds")
graphics::abline(h = capacity / 1e6, col = "#d7301f", lty = 2, lwd = 2)
graphics::abline(h = cap_naive / 1e6, col = "grey40", lty = 3, lwd = 2)
legend("topleft", bty = "n", lwd = c(2.5, 2, 2), lty = c(1, 2, 3),
       col = c("#238b45", "#d7301f", "grey40"),
       legend = c("living population", "calc_capacity_cdr (allocated, with squash)",
                  "calc_capacity (needed without reclaim)"))
close_png()

# ── plot 2: active slots (squash sawtooth) vs living population vs capacity ──────────
open_png(file.path(out_dir, "long_run_squash_utilization.png"), width = 1050, height = 680, res = 110)
matplot(yrs, cbind(active = active / 1e6, living = live_pop / 1e6), type = "l", lty = 1, lwd = 2,
        col = c("#2c7fb8", "#238b45"), ylim = c(0, capacity / 1e6 * 1.05),
        xlab = "year", ylab = "agents (millions)",
        main = "Active slots reclaimed by annual squash() stay near the living count")
graphics::abline(h = capacity / 1e6, col = "#d7301f", lty = 2, lwd = 2)
legend("topleft", bty = "n", lwd = c(2, 2, 2), lty = c(1, 1, 2),
       col = c("#2c7fb8", "#238b45", "#d7301f"),
       legend = c("active slots (people$count; resets each year at squash)",
                  "living population", "allocated capacity"))
close_png()

# ── plot 3: daily births and deaths ─────────────────────────────────────────────────
open_png(file.path(out_dir, "long_run_squash_vitals.png"), width = 1050, height = 680, res = 110)
matplot(yrs[-1], cbind(births = births_d, deaths = deaths_d), type = "l", lty = 1, lwd = 1.5,
        col = c("#2c7fb8", "#d7301f"), xlab = "year", ylab = "events / day",
        main = sprintf("Daily births (CBR %d) and deaths (CDR %d) over the growing population", cbr, cdr))
legend("topleft", legend = c("births", "deaths"), bty = "n", lwd = 1.5,
       col = c("#2c7fb8", "#d7301f"))
close_png()

if (to_png) cat(sprintf("Plots written to: %s\n", normalizePath(out_dir)))
