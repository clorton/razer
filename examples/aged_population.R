#!/usr/bin/env Rscript

# End-to-end demographics example: build a realistic age-structured population and give
# every agent a realistic date of death. Two razer demographics pieces (ported from
# laser.core) do the work:
#   * aliased_distribution() (via the sample_pyramid_ages() helper) — draw each agent's
#     age from a population pyramid, weighted by the population in each age band;
#   * KaplanMeierEstimator — sample each agent's age AT DEATH from a life table,
#     conditioned on the age they are alive at now (so nobody dies in the past).
#
# Run from anywhere:  Rscript examples/aged_population.R

# `library(pkg)` attaches a package so its exported names resolve unqualified.
library(razer)

# Reproducible R-side draws (the within-band age jitter and the histogram). NOTE: the
# alias bin-selection and the Kaplan–Meier sampling use the package's internal,
# thread-local RNG, which is NOT affected by set.seed — so those parts vary run to run.
set.seed(20240605L)

n_agents <- 1000000L   # one million agents

# ── 1. realistic ages from a population pyramid ─────────────────────────────────────
# Resolve this script's directory so the data file and output plot resolve no matter
# the working directory. `commandArgs(FALSE)` includes Rscript's own `--file=<path>`.
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else "."

# `load_pyramid_csv` parses the "Age,M,F" pyramid into a start/end/M/F integer matrix.
pyramid <- load_pyramid_csv(file.path(script_dir, "data", "pyramid_example.csv"))

# `sample_pyramid_ages` builds an aliased_distribution over the per-band populations
# (M + F), draws a band per agent in proportion to its population, then a uniform day
# within the band's year range — returning ages in DAYS.
ages_days  <- sample_pyramid_ages(pyramid, n_agents)
ages_years <- ages_days / 365

# ── 2. a synthetic life table (cumulative deaths by year) ───────────────────────────
# The Kaplan–Meier estimator is built from CUMULATIVE deaths by year for a synthetic
# cohort. We synthesize one from a per-year mortality hazard q(age): a small baseline,
# an exponential (Gompertz-like) rise with age, and an elevated infant rate. Then we run
# a cohort of 100,000 newborns through it, recording deaths each year, and force any
# survivors to die in the final year so the table is complete. The slope is chosen so
# the hazard approaches 1 by ~age 100, leaving essentially no survivors to pile up in
# the forced final year (implied life expectancy at birth ~73 years).
n_years <- 106L                                    # ages 0..105 (max_year = 105)
age     <- 0:(n_years - 1L)
q       <- 0.0004 + 1e-5 * exp(0.115 * age)        # Gompertz-Makeham-ish hazard
q[1]    <- 0.02                                     # infant (age 0) mortality
q       <- pmin(q, 1)

cohort <- 100000
deaths <- numeric(n_years)
alive  <- cohort
for (a in seq_len(n_years)) {
  # `if (...) a else b` is an expression; the final year takes everyone still alive.
  d         <- if (a == n_years) alive else round(alive * q[a])
  deaths[a] <- d
  alive     <- alive - d
}
cumulative_deaths <- cumsum(deaths)                # non-decreasing, ends at `cohort`

# Build the estimator. `kaplan_meier_estimator` validates the table and prepends its
# own internal leading zero.
km <- kaplan_meier_estimator(cumulative_deaths)

# ── 3. assign each agent a date (age) of death ──────────────────────────────────────
# `predict_age_at_death` samples an age at death IN DAYS for each agent, conditioned on
# its current age (the result is never earlier than the agent's current age). max_year
# = -1L uses the last year in the life table (100).
timing <- system.time(
  death_age_days <- km$predict_age_at_death(ages_days, -1L)
)
death_age_years    <- death_age_days / 365
remaining_lifespan <- (death_age_days - ages_days) / 365   # years of life left

# Sanity checks + summary. `stopifnot` errors if any condition is FALSE.
stopifnot(all(death_age_days >= ages_days))                # nobody dies in the past
stopifnot(all(death_age_days < n_years * 365))             # within the life table
cat(sprintf("aged_population: %s agents; ages %.1f-%.1f yr (mean %.1f).\n",
            format(n_agents, big.mark = ","),
            min(ages_years), max(ages_years), mean(ages_years)))
cat(sprintf("  predicted age at death: mean %.1f yr (range %.1f-%.1f); ",
            mean(death_age_years), min(death_age_years), max(death_age_years)))
cat(sprintf("mean remaining lifespan %.1f yr.\n", mean(remaining_lifespan)))
# Life expectancy at birth implied by the table = mean age at death of newborns.
e0 <- mean(km$predict_age_at_death(integer(100000L), -1L)) / 365
cat(sprintf("  implied life expectancy at birth e0 = %.1f yr; prediction took %.3f s.\n",
            e0, timing[["elapsed"]]))

# ── 4. plot the population age structure and the age-at-death distribution ──────────
# Device-aware: write a PNG when run non-interactively (Rscript); draw to the active
# device (e.g. RStudio's Plots pane) when sourced interactively.
to_png    <- !interactive()
open_png  <- function(path, ...) if (to_png) grDevices::png(path, ...)
close_png <- function() if (to_png) grDevices::dev.off()
out_dir <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
plot_path <- file.path(out_dir, "aged_population.png")
open_png(plot_path, width = 1100, height = 850, res = 120)
op <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4.5, 2.5, 1))

# Top: the generated age distribution (the pyramid we sampled into being).
breaks <- seq(0, n_years, by = 1)
hist(ages_years, breaks = breaks, col = "steelblue", border = "white",
     xlab = "age (years)", ylab = "agents",
     main = sprintf("Generated population age structure (%s agents)",
                    format(n_agents, big.mark = ",")))

# Bottom: the predicted age-at-death distribution — note the infant-mortality spike at
# the left and the bulk of deaths concentrated in old age.
hist(death_age_years, breaks = breaks, col = "firebrick", border = "white",
     xlab = "age at death (years)", ylab = "agents",
     main = "Predicted age at death (Kaplan-Meier, conditioned on current age)")
graphics::abline(v = e0, col = "grey30", lty = 2)            # life expectancy at birth
graphics::legend("topleft", legend = sprintf("e0 = %.1f yr", e0), bty = "n", lty = 2)

graphics::par(op)
close_png()
if (to_png) cat(sprintf("wrote age / age-at-death plot to %s\n", plot_path))
