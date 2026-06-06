# Tests for KaplanMeierEstimator: given a cumulative-deaths-by-year life table, sample
# a year (or age in days) of death conditioned on survival to a given age. The RNG is
# thread-local (not R-seedable), so most checks are bounds/monotonicity/statistical.
# A life table that concentrates all deaths in a single year gives a deterministic
# year-of-death, used for an exact check. Written given-when-then.

# All 100 deaths happen during year 50: cumulative is 0 through year 49, then 100
# (this IS the cumulative table — 0 for years 0..49, then a flat 100 for years 50..100).
deaths_in_year_50 <- c(rep(0, 50), rep(100, 51))           # length 101 (years 0..100)
# A smoothly increasing table (10 deaths/yr for 80 yr, then 100/yr) for distribution checks.
smooth_table      <- cumsum(c(rep(10, 80), rep(100, 21)))  # length 101

test_that("predict_year_of_death is exact when all deaths fall in one year", {
  # Given a table where the entire cohort dies during year 50
  # When predicting the year of death for individuals aged 0, 25, and 50
  # Then every prediction is exactly year 50 (the only year with deaths), regardless
  #      of current age. Failure would mean the conditional draw or the binary-search
  #      year lookup is off by one or mis-conditioned.
  km <- kaplan_meier_estimator(deaths_in_year_50)

  expect_true(all(km$predict_year_of_death(rep(0L,  1000L), -1L) == 50L))
  expect_true(all(km$predict_year_of_death(rep(25L, 1000L), -1L) == 50L))
  expect_true(all(km$predict_year_of_death(rep(50L, 1000L), -1L) == 50L))
})

test_that("predict_year_of_death stays in [age, max_year] and conditions on age", {
  # Given a smooth life table
  # When predicting from a young age vs. an old age
  # Then all predictions are >= the current age and <= max_year, and the mean year of
  #      death is LATER for older survivors (survival conditioning). Failure would mean
  #      deaths are sampled unconditionally or escape the [age, max_year] range.
  km <- kaplan_meier_estimator(smooth_table)              # max_year = 100

  young <- km$predict_year_of_death(rep(10L, 50000L), -1L)
  old   <- km$predict_year_of_death(rep(70L, 50000L), -1L)

  expect_true(all(young >= 10L & young <= 100L))
  expect_true(all(old   >= 70L & old   <= 100L))
  expect_gt(mean(old), mean(young))                       # older -> later death
})

test_that("predict_year_of_death respects an explicit max_year cap", {
  # Given a smooth table and an explicit max_year of 60
  # When predicting from age 40
  # Then no predicted year exceeds 60. Failure would mean the cap is ignored.
  km <- kaplan_meier_estimator(smooth_table)
  yod <- km$predict_year_of_death(rep(40L, 20000L), 60L)
  expect_true(all(yod >= 40L & yod <= 60L))
})

test_that("predict_age_at_death returns days in [current age, (max_year+1)*365)", {
  # Given a smooth table and individuals aged 40 years (in days)
  # When predicting the age at death in days
  # Then every result is >= the current age (death is never in the past) and strictly
  #      below the (max_year+1)*365 ceiling. With the single-year-50 table the year of
  #      death is fixed, so the day of death lands in [50*365, 51*365). Failure would
  #      mean the day-within-year logic or the year offset is wrong.
  km <- kaplan_meier_estimator(smooth_table)
  ages_days <- rep(40L * 365L, 50000L)
  aad <- km$predict_age_at_death(ages_days, -1L)
  expect_true(all(aad >= ages_days))                      # not earlier than now
  expect_true(all(aad <  101L * 365L))                    # below (100+1)*365

  km50 <- kaplan_meier_estimator(deaths_in_year_50)
  aad50 <- km50$predict_age_at_death(rep(0L, 20000L), -1L)
  expect_true(all(aad50 >= 50L * 365L & aad50 < 51L * 365L))
})

test_that("cumulative_deaths round-trips the source (without the internal leading 0)", {
  # Given a source table
  # When read back via cumulative_deaths()
  # Then it equals the input (the internally-prepended zero is excluded). Failure would
  #      mean the accessor leaks the internal padding or mis-stores the table.
  km <- kaplan_meier_estimator(smooth_table)
  expect_equal(km$cumulative_deaths(), as.numeric(smooth_table))
})

test_that("kaplan_meier_estimator and predictions reject invalid input", {
  # Given non-monotonic / negative tables and out-of-range query ages
  # When the constructor or a predictor is called
  # Then it errors rather than producing nonsense.
  # Failure would risk out-of-bounds draws or silently wrong lifespans.
  expect_error(kaplan_meier_estimator(c(10, 5, 20)), "non-decreasing")
  expect_error(kaplan_meier_estimator(c(-1, 2, 3)), "non-negative")

  km <- kaplan_meier_estimator(smooth_table)              # 101 years (max_year 100)
  expect_error(km$predict_year_of_death(c(101L), -1L), "max_year")        # age past table
  expect_error(km$predict_age_at_death(c(101L * 365L), -1L), "ages in days") # day past ceiling
})
