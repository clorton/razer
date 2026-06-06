#!/usr/bin/env Rscript

# A century-long demographic run that stays within a bounded agent array by reclaiming the
# slots of the dead with squash(). This is the use case calc_capacity_cdr() exists for.
#
# Setup: one node, 1,000,000 initial agents, run for 100 years (daily steps) under a crude
# BIRTH rate of 30 / 1,000 / year and a crude DEATH rate of 15 / 1,000 / year (so the
# population grows at ~15 / 1,000 / year net). With no reclaim you would have to allocate a
# slot for every agent EVER born — tens of millions over a century (what calc_capacity()
# estimates). Instead we:
#
#   * size the array with calc_capacity_cdr(CBR, CDR), which bounds the PEAK LIVING
#     population (net births minus an underestimated death rate) rather than cumulative
#     births, and
#   * call squash() once a year to compact the living agents to the front of every per-agent
#     Column and free the dead slots, so the active count stays near the living count and
#     well under the allocated capacity.
#
# Births enter the maternal-immunity compartment M (the births() kernel draws a maternal
# timer and a Kaplan-Meier date of death per newborn); step_si wanes M -> S. Mortality is
# age-independent ("crude"): every agent's lifespan is exponential with the CDR hazard, via
# a Kaplan-Meier life table built from an exponential survival curve, so the realized death
# rate is ~CDR regardless of age structure.
#
# Run from anywhere:  Rscript examples/long_run_squash.R   (takes a few minutes)
# Output PNGs are written next to this script in examples/output/.

library(razer)

set.seed(20260606L)
states <- laser_states()                 # c(S = 0, E = 1, I = 2, R = 3, M = 4, D = -1)

# ── parameters ────────────────────────────────────────────────────────────────────
N0        <- 1000000L                     # initial living population
years     <- 100L
nticks    <- years * 365L                 # daily steps
cbr       <- 30                           # crude birth rate (per 1,000 per year)
cdr       <- 15                           # crude death rate (per 1,000 per year)
maternal_duration <- dist_normal(270, 20) # M -> S waning (~9 months)
squash_every <- 365L                      # reclaim dead slots once a year

# ── capacity: peak-living bound (with squash) vs cumulative-births bound (without) ──
# Both rate grids are constant over the run. calc_capacity_cdr is what we ALLOCATE;
# calc_capacity is shown only to quantify how many slots squashing saves.
br_grid   <- matrix(cbr, nrow = nticks - 1L, ncol = 1L)
dr_grid   <- matrix(cdr, nrow = nticks - 1L, ncol = 1L)
capacity  <- as.integer(calc_capacity_cdr(br_grid, dr_grid, N0, safety_factor = 1))
cap_naive <- as.integer(min(calc_capacity(br_grid, N0, safety_factor = 1), .Machine$integer.max))
cat(sprintf("capacity (squash-aware calc_capacity_cdr) = %s slots\n", format(capacity, big.mark = ",")))
cat(sprintf("capacity (no-reclaim calc_capacity)       = %s slots  (%.1fx more)\n",
            format(cap_naive, big.mark = ","), cap_naive / capacity))

# ── exponential life table for crude (age-independent) mortality at the CDR ─────────
# Constant hazard h = -log(1 - CDR/1000) per year gives an exponential lifespan, so the
# per-capita death rate is ~CDR at every age. Build the cumulative-deaths curve a KM
# estimator samples ages at death from (used for both initial agents and newborns).
max_age   <- 200L
h_year    <- -log(1 - cdr / 1000)
surv      <- exp(-h_year * (0:max_age))           # l(a) = exp(-h a)
cohort    <- 1e7
cumulative_deaths <- round((1 - surv) * cohort)   # life table: deaths by age (memoryless)
km        <- kaplan_meier_estimator(cumulative_deaths)

# ── people: per-agent arrays allocated to the squash-aware capacity ─────────────────
people <- new.env()
people$count    <- N0
people$capacity <- capacity
people$state    <- allocate_scalar("u8",  capacity)
people$nodeid   <- allocate_scalar("u16", capacity)
people$timer    <- allocate_scalar("u16", capacity)
people$dob      <- allocate_scalar("i32", capacity)   # date of birth (negative = born before t0)
people$dod      <- allocate_scalar("u32", capacity)   # absolute tick of death

# Initial agents: ages from the exponential stationary structure (younger more common),
# everyone starts Susceptible; reserved slots (state 0 = S, nodeid 0) are inert until a
# birth activates them. dob = -age; dod = current tick (0) + a memoryless remaining lifespan
# drawn (conditioned on the current age) from the KM life table.
# Cap initial ages to 90 years: the exponential tail can exceed the life table's range,
# and centenarians are vanishingly rare anyway.
age_days     <- pmin(as.integer(floor(stats::rexp(N0, rate = h_year / 365))), 90L * 365L)
age_at_death <- km$predict_age_at_death(age_days, -1L)   # -1L: full life table
people$state$set(rep(states[["S"]], capacity))
people$nodeid$set(rep(0L, capacity))
people$dob$set(c(-age_days, rep(0L, capacity - N0)))
people$dod$set(c(as.integer(age_at_death - age_days), rep(0L, capacity - N0)))

# ── nodes: a single-node census (M and S; no disease) + birth/death flows ───────────
n_nodes <- 1L
nodes <- new.env()
nodes$M <- allocate_vector("i32", nticks, n_nodes)
nodes$S <- allocate_vector("i32", nticks, n_nodes)
nodes$N <- allocate_vector("i32", nticks, n_nodes)
zeros <- rep(0L, (nticks - 1L) * n_nodes)
nodes$M$set(c(0L,  zeros))
nodes$S$set(c(N0,  zeros))
nodes$N$set(c(N0,  zeros))
nodes$birth_rate <- values_map(cbr / 1000 / 365, nticks, n_nodes)   # per person per day

# Per-tick series we record for the plots.
live_pop  <- numeric(nticks)        # living agents (M + S census)  -- one per recorded tick
active    <- numeric(nticks)        # people$count (live + not-yet-squashed dead slots)
births_d  <- numeric(nticks - 1L)
deaths_d  <- numeric(nticks - 1L)
live_pop[1] <- N0; active[1] <- N0
placeholder <- dist_constant(7)     # step_si needs an inf-duration arg; unused (no E here)

# ── the 100-year loop ───────────────────────────────────────────────────────────────
pb <- utils::txtProgressBar(min = 0L, max = nticks - 1L, style = 3)
every <- max(1L, (nticks - 1L) %/% 100L)
for (tick in seq_len(nticks - 1L)) {
  t <- tick - 1L
  carry_forward_states(list(nodes$M, nodes$S), t, total = nodes$N)
  # M -> S waning (step_si also does E -> I, but there is no E here; onset is 0).
  prog <- step_si(people$state, people$timer, people$nodeid, people$count, n_nodes, placeholder)
  move_count(nodes$M, nodes$S, prog$waned, t)
  # Births: CBR newborns into M, each with a maternal timer and a KM date of death.
  b <- births(people$state, people$timer, people$nodeid, people$dob, people$dod,
              people$count, n_nodes, nodes$birth_rate, maternal_duration, km, t)
  people$count <- b$count
  move_count(NULL, nodes$M, b$born, t)
  births_d[tick] <- sum(b$born)
  # Crude mortality: retire agents whose dod has arrived; decrement their compartments.
  d <- mortality(people$state, people$dod, people$nodeid, people$count, n_nodes, t)
  move_count(nodes$M, NULL, d$m, t); move_count(nodes$S, NULL, d$s, t)
  deaths_d[tick] <- sum(d$m + d$s)
  # Once a year, reclaim the slots of the dead so the active count tracks the LIVING count
  # instead of growing toward the cumulative-births total. squash() compacts every
  # per-agent Column by the same alive mask and updates people$count; the census is
  # aggregate and untouched.
  if (tick %% squash_every == 0L) people$count <- squash(people)
  live_pop[tick + 1L] <- sum(nodes$M$values()[tick + 1L, ], nodes$S$values()[tick + 1L, ])
  active[tick + 1L]   <- people$count
  if (tick %% every == 0L || tick == nticks - 1L) utils::setTxtProgressBar(pb, tick)
}
close(pb)

# ── report ────────────────────────────────────────────────────────────────────────
person_years <- sum(live_pop) / 365
cat(sprintf("\nliving population: %s -> %s over %d years (net growth %.2fx)\n",
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
out_dir    <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

yrs <- (seq_len(nticks) - 1L) / 365

# ── plot 1: living population vs the two capacity bounds ─────────────────────────────
grDevices::png(file.path(out_dir, "long_run_squash_population.png"), width = 1050, height = 680, res = 110)
plot(yrs, live_pop / 1e6, type = "l", lwd = 2.5, col = "#238b45",
     ylim = c(0, cap_naive / 1e6 * 1.05), xlab = "year", ylab = "agents (millions)",
     main = "100-year run: living population vs. agent-array capacity bounds")
graphics::abline(h = capacity / 1e6, col = "#d7301f", lty = 2, lwd = 2)
graphics::abline(h = cap_naive / 1e6, col = "grey40", lty = 3, lwd = 2)
legend("topleft", bty = "n", lwd = c(2.5, 2, 2), lty = c(1, 2, 3),
       col = c("#238b45", "#d7301f", "grey40"),
       legend = c("living population", "calc_capacity_cdr (allocated, with squash)",
                  "calc_capacity (needed without reclaim)"))
grDevices::dev.off()

# ── plot 2: active slots (squash sawtooth) vs living population vs capacity ──────────
grDevices::png(file.path(out_dir, "long_run_squash_utilization.png"), width = 1050, height = 680, res = 110)
matplot(yrs, cbind(active = active / 1e6, living = live_pop / 1e6), type = "l", lty = 1, lwd = 2,
        col = c("#2c7fb8", "#238b45"), ylim = c(0, capacity / 1e6 * 1.05),
        xlab = "year", ylab = "agents (millions)",
        main = "Active slots reclaimed by annual squash() stay near the living count")
graphics::abline(h = capacity / 1e6, col = "#d7301f", lty = 2, lwd = 2)
legend("topleft", bty = "n", lwd = c(2, 2, 2), lty = c(1, 1, 2),
       col = c("#2c7fb8", "#238b45", "#d7301f"),
       legend = c("active slots (people$count; resets each year at squash)",
                  "living population", "allocated capacity"))
grDevices::dev.off()

# ── plot 3: daily births and deaths ─────────────────────────────────────────────────
grDevices::png(file.path(out_dir, "long_run_squash_vitals.png"), width = 1050, height = 680, res = 110)
matplot(yrs[-1], cbind(births = births_d, deaths = deaths_d), type = "l", lty = 1, lwd = 1.5,
        col = c("#2c7fb8", "#d7301f"), xlab = "year", ylab = "events / day",
        main = sprintf("Daily births (CBR %d) and deaths (CDR %d) over the growing population", cbr, cdr))
legend("topleft", legend = c("births", "deaths"), bty = "n", lwd = 1.5,
       col = c("#2c7fb8", "#d7301f"))
grDevices::dev.off()

cat(sprintf("Plots written to: %s\n", normalizePath(out_dir)))
