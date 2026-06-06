// ════════════════════════════════════════════════════════════════════════════
// kmestimator.rs — Kaplan–Meier sampling of a (remaining) lifespan.
//
// `KaplanMeierEstimator` turns a CUMULATIVE-DEATHS-by-year life table into a sampler
// for "given an individual alive at age X, when do they die?". It is the engine for
// assigning each agent a realistic date of death at creation, conditioned on the age
// they are created at. A faithful port of laser.core's `KaplanMeierEstimator`.
//
// Input: `cumulative_deaths[y]` = number of a synthetic cohort dead by the end of year
// `y` (non-decreasing). We prepend a 0 so index `y` means "deaths strictly before year
// y". To draw a year of death for someone currently at age `a` (and capped at
// `max_year`), we pick a death-RANK uniformly among the cohort members not yet dead at
// age `a` — i.e. in `(cumulative[a], cumulative[max_year+1]]` — and read off which year
// that rank falls in (a binary search). Conditioning on survival to age `a` is exactly
// the Kaplan–Meier conditional-survival idea.
//
// Orientation for readers coming from C / C++ / C# / Python:
//   * `#[extendr]` exposes the struct as an opaque R handle and the `impl` methods as
//     `obj$method(...)`; the free `kaplan_meier_estimator(...)` fn is the constructor.
//   * `slice.partition_point(|&x| x < draw)` is C++'s `lower_bound`: on a sorted slice
//     it returns the index of the first element `>= draw` (NumPy `searchsorted` left).
//   * `gen_range(lo..=hi)` is an INCLUSIVE uniform integer draw; `lo..hi` is half-open.
//   * Work is split across cores with `par_chunks_mut` (Rayon) + a per-chunk RNG from
//     `rng.rs` (the same across-agents parallelism as the simulation kernels, and like
//     the `@njit(parallel=True)` original). The RNG is seedable via `set_seed`.
//   * A negative `max_year` argument is the sentinel for "use the last year in the
//     table" (extendr methods have no native default-argument support).
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;
use rayon::prelude::*;
use rand::Rng;
use crate::rng;

const DAYS_PER_YEAR: i64 = 365;

/// A Kaplan–Meier sampler over a cumulative-deaths-by-year life table.
///
/// Construct with [kaplan_meier_estimator()] from a non-decreasing vector of
/// cumulative deaths by year. The handle is opaque to R.
///
/// @export
#[extendr]
pub struct KaplanMeierEstimator {
    // Cumulative deaths with a leading 0 prepended: `cd[y]` = deaths strictly before
    // year `y`, for y in 0..=n_years. Length is `n_years + 1`. i64 holds national-
    // scale cohort sizes and keeps the binary-search arithmetic exact.
    cd: Vec<i64>,
}

impl KaplanMeierEstimator {
    fn n_years(&self) -> usize {
        self.cd.len() - 1
    }

    // Resolve the max-year argument: a negative value means "last year in the table".
    fn resolve_max_year(&self, max_year: i32) -> usize {
        let n_years = self.n_years();
        if max_year < 0 {
            n_years - 1
        } else {
            let m = max_year as usize;
            assert!(
                m < n_years,
                "max_year ({m}) must be less than the number of years in the table ({n_years})"
            );
            m
        }
    }

    // Draw a year of death for an individual currently at integer age `age_years`,
    // capped at `max_year`. `total_deaths` is `cd[max_year + 1]` (precomputed once).
    fn draw_year_of_death<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
        age_years: usize,
        max_year: usize,
        total_deaths: i64,
    ) -> usize {
        let already_deceased = self.cd[age_years];
        if already_deceased >= total_deaths {
            // The whole cohort is gone by `max_year`; die in the final year.
            return max_year;
        }
        // Uniform death-rank among those still alive at `age_years`, then locate its
        // year: first index whose cumulative count reaches `draw`, minus one.
        let draw = rng.gen_range(already_deceased + 1..=total_deaths);
        self.cd.partition_point(|&x| x < draw) - 1
    }

    // Sample one age AT death (in days) for an individual currently `age_days` old:
    // draw a year of death, then a day within it (a uniform day of a later year, or — if
    // death is this same year — a day at or after the current day-of-year). `pub(crate)`
    // so kernels (e.g. births) can assign a newborn's lifespan without going through R.
    pub(crate) fn sample_age_at_death_days<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
        age_days: i64,
        max_year: usize,
        total_deaths: i64,
    ) -> i64 {
        let age_years = (age_days / DAYS_PER_YEAR) as usize;
        let yod = self.draw_year_of_death(rng, age_years, max_year, total_deaths);
        let doy = if (age_days / DAYS_PER_YEAR) < yod as i64 {
            rng.gen_range(0..DAYS_PER_YEAR)          // dies in a later year: any day
        } else {
            let age_doy = age_days % DAYS_PER_YEAR;  // dies this year: today or later
            rng.gen_range(age_doy..DAYS_PER_YEAR)
        };
        yod as i64 * DAYS_PER_YEAR + doy
    }

    // Convenience for a NEWBORN (age 0), using the full life table. Returns the age at
    // death in days (equivalently the days until death for someone born now).
    pub(crate) fn sample_newborn_age_at_death<R: Rng + ?Sized>(&self, rng: &mut R) -> i64 {
        let max_year = self.n_years() - 1;
        let total_deaths = self.cd[max_year + 1];
        self.sample_age_at_death_days(rng, 0, max_year, total_deaths)
    }

    // Per-chunk parallel helper shared by the two predict methods. Seeds a per-chunk RNG
    // from the shared seed (fixed chunk size for reproducibility); `f` receives the chunk,
    // its starting global index, and that RNG.
    fn par_fill(&self, len: usize, f: impl Fn(&mut [i32], usize, &mut rng::ModelRng) + Sync) -> Vec<i32> {
        let mut out = vec![0i32; len];
        let base = rng::next_call_base();
        let chunk = rng::RNG_CHUNK;
        out.par_chunks_mut(chunk)
            .enumerate()
            .for_each(|(ci, c)| {
                let mut r = rng::chunk_rng(base, ci);
                f(c, ci * chunk, &mut r);
            });
        out
    }
}

#[extendr]
impl KaplanMeierEstimator {
    /// Predict a year of death for each individual given their current age in years.
    ///
    /// For each age, samples a year of death `>= age` and `<= max_year`, conditioned on
    /// survival to that age (Kaplan–Meier). Ages must be in `0..=max_year`.
    ///
    /// @param ages_years Integer vector of current ages in whole years.
    /// @param max_year   Maximum year of death to consider; pass a negative value
    ///   (e.g. `-1L`) to use the last year in the life table.
    /// @return An integer vector of predicted years of death (same length as input).
    /// @examples
    /// km <- kaplan_meier_estimator(cumsum(c(rep(10, 80), rep(100, 21))))
    /// km$predict_year_of_death(c(40L, 50L, 60L), -1L)
    /// @export
    fn predict_year_of_death(&self, ages_years: Vec<i32>, max_year: i32) -> Vec<i32> {
        let max_year = self.resolve_max_year(max_year);
        let total_deaths = self.cd[max_year + 1];
        for &a in &ages_years {
            assert!(a >= 0, "ages must be non-negative, got {a}");
            assert!(
                (a as usize) <= max_year,
                "all ages must be <= max_year ({max_year}), got {a}"
            );
        }
        self.par_fill(ages_years.len(), |c, start, rng| {
            for (j, slot) in c.iter_mut().enumerate() {
                let age = ages_years[start + j] as usize;
                *slot = self.draw_year_of_death(rng, age, max_year, total_deaths) as i32;
            }
        })
    }

    /// Predict an age at death (in DAYS) for each individual given their age in days.
    ///
    /// Samples the year of death as in [predict_year_of_death()], then a day within
    /// that year: a uniform day of a later year, or — if death falls in the individual's
    /// current year — a uniform day at or after their current day-of-year (so the
    /// predicted age at death is never earlier than the current age). Ages in days must
    /// be `< (max_year + 1) * 365`.
    ///
    /// @param ages_days Integer vector of current ages in whole days.
    /// @param max_year  Maximum year of death to consider; pass a negative value
    ///   (e.g. `-1L`) to use the last year in the life table.
    /// @return An integer vector of predicted ages at death in days (same length).
    /// @examples
    /// km <- kaplan_meier_estimator(cumsum(c(rep(10, 80), rep(100, 21))))
    /// km$predict_age_at_death(c(40L, 50L, 60L) * 365L, -1L)
    /// @export
    fn predict_age_at_death(&self, ages_days: Vec<i32>, max_year: i32) -> Vec<i32> {
        let max_year = self.resolve_max_year(max_year);
        let total_deaths = self.cd[max_year + 1];
        let limit = (max_year as i64 + 1) * DAYS_PER_YEAR;
        for &ad in &ages_days {
            assert!(ad >= 0, "ages must be non-negative, got {ad}");
            assert!(
                (ad as i64) < limit,
                "all ages in days must be < (max_year + 1) * 365 = {limit}, got {ad}"
            );
        }
        self.par_fill(ages_days.len(), |c, start, rng| {
            for (j, slot) in c.iter_mut().enumerate() {
                let age_days = ages_days[start + j] as i64;
                *slot = self.sample_age_at_death_days(rng, age_days, max_year, total_deaths) as i32;
            }
        })
    }

    /// The cumulative-deaths-by-year table (without the internal leading zero).
    /// @return A numeric vector of length equal to the number of years.
    /// @export
    fn cumulative_deaths(&self) -> Vec<f64> {
        self.cd[1..].iter().map(|&x| x as f64).collect()
    }
}

/// Build a [KaplanMeierEstimator] from cumulative deaths by year.
///
/// `cumulative_deaths[y]` is the number of a synthetic cohort dead by the end of year
/// `y`; it must be non-negative and monotonically non-decreasing. Values are rounded to
/// whole numbers. (A leading zero is prepended internally; do not include it yourself.)
///
/// @param cumulative_deaths A non-decreasing numeric vector of cumulative deaths by
///   year (length >= 1).
/// @return A `KaplanMeierEstimator` object.
/// @examples
/// # toy life table: 10 deaths/year for 80 years, then 100/year for 21 more
/// km <- kaplan_meier_estimator(cumsum(c(rep(10, 80), rep(100, 21))))
/// @export
#[extendr]
fn kaplan_meier_estimator(cumulative_deaths: Vec<f64>) -> KaplanMeierEstimator {
    let n = cumulative_deaths.len();
    assert!(n >= 1, "`cumulative_deaths` must have at least one year");

    let src: Vec<i64> = cumulative_deaths
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            assert!(
                c.is_finite() && c >= 0.0,
                "cumulative_deaths[{i}] must be finite and non-negative, got {c}"
            );
            c.round() as i64
        })
        .collect();

    for w in src.windows(2) {
        assert!(
            w[1] >= w[0],
            "`cumulative_deaths` must be monotonically non-decreasing"
        );
    }

    // Prepend the leading zero: cd[0] = 0, cd[y+1] = cumulative_deaths[y].
    let mut cd = Vec::with_capacity(n + 1);
    cd.push(0);
    cd.extend_from_slice(&src);
    KaplanMeierEstimator { cd }
}

extendr_module! {
    mod kmestimator;
    impl KaplanMeierEstimator;
    fn kaplan_meier_estimator;
}
