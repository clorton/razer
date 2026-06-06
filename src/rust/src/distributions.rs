use extendr_api::prelude::*;
use rand::Rng;
use rand::distributions::Distribution as RandDistribution;
use rand::distributions::{Open01, Uniform};
use rand_distr::{Beta, Exp, Gamma, LogNormal, Normal, Poisson};

// ── Distribution families ───────────────────────────────────────────────────────
//
// `DistKind` is the closed set of supported families.  Each variant holds only a
// small, immutable parameter struct (`Send + Sync`), so a `Distribution` can be
// shared by reference across Rayon worker threads and sampled concurrently with no
// locking and no allocation. Sampling matches on `&self.kind` and borrows the
// inner sampler, so no variant is required to be `Copy`.

#[derive(Clone)]
enum DistKind {
    /// Degenerate distribution: every draw returns the same value.
    Constant(f64),
    /// Normal (Gaussian) distribution, stored as (mean, std_dev).
    Normal(Normal<f64>),
    /// Continuous uniform distribution on a half-open interval [low, high).
    Uniform(Uniform<f64>),
    /// Gamma distribution, stored as (shape k, scale θ).
    Gamma(Gamma<f64>),
    /// Poisson distribution with rate λ; draws are non-negative integer counts.
    Poisson(Poisson<f64>),
    /// Beta distribution on (0, 1), stored as (α, β).
    Beta(Beta<f64>),
    /// Exponential distribution with rate λ (mean 1/λ).
    Exp(Exp<f64>),
    /// Logistic distribution with location μ and scale s. `rand_distr` has no
    /// logistic sampler, so it is drawn by inverse-CDF transform of a uniform.
    Logistic { location: f64, scale: f64 },
    /// Log-normal distribution whose log is Normal(meanlog, sdlog).
    LogNormal(LogNormal<f64>),
}

/// A parameterized probability distribution that can be sampled repeatedly.
///
/// Build one with a family constructor such as `dist_normal()` or `dist_constant()`,
/// then either sample it from Rust via `sample()` (Pattern B: the caller owns the
/// RNG and passes it in), or from R via `$sample_one()` / `$sample_n()`.
///
/// Draws are always floating-point. Callers that need integer values (e.g. a
/// whole-tick state timer) are responsible for rounding or truncating as
/// appropriate — the simulation kernels round to the nearest tick.
///
/// The handle is opaque to R — it is passed to simulation kernels (e.g.
/// `step_sir` or `transmission`) by reference, so the same object can be reused
/// every tick and shared across all worker threads.
///
/// @export
#[extendr]
pub struct Distribution {
    kind: DistKind,
}

impl Distribution {
    /// Draw a single `f64` sample using the supplied RNG.
    ///
    /// This is the Pattern B entry point used by simulation kernels: the caller
    /// owns the RNG (typically a per-thread `thread_rng`) and passes a mutable
    /// reference in, so one shared `&Distribution` can be sampled concurrently
    /// from every Rayon worker.
    pub fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        match &self.kind {
            DistKind::Constant(v) => *v,
            DistKind::Normal(d) => d.sample(rng),
            DistKind::Uniform(d) => d.sample(rng),
            DistKind::Gamma(d) => d.sample(rng),
            DistKind::Poisson(d) => d.sample(rng),
            DistKind::Beta(d) => d.sample(rng),
            DistKind::Exp(d) => d.sample(rng),
            DistKind::Logistic { location, scale } => {
                // Inverse-CDF transform: X = μ + s·logit(U), U ~ Uniform(0, 1).
                // Open01 excludes 0 and 1, so logit(U) is always finite.
                let u: f64 = rng.sample(Open01);
                *location + *scale * (u / (1.0 - u)).ln()
            }
            DistKind::LogNormal(d) => d.sample(rng),
        }
    }
}

#[extendr]
impl Distribution {
    /// Draw a single sample using a thread-local RNG.
    ///
    /// Convenience for interactive use. Simulation kernels do not call this — they
    /// use the internal sampler with an explicit, reusable RNG for performance.
    ///
    /// @return A single numeric (double) draw from the distribution.
    /// @export
    fn sample_one(&self) -> f64 {
        self.sample(&mut crate::rng::single_rng())
    }

    /// Draw `n` samples using a thread-local RNG, returned as a numeric vector.
    ///
    /// Drawing a whole batch in one call avoids per-sample R↔Rust overhead, which
    /// makes it practical to validate the sampler against large empirical samples
    /// (e.g. one million draws) from R.
    ///
    /// @param n Number of samples to draw; must be non-negative.
    /// @return A numeric (double) vector of length `n`.
    /// @export
    fn sample_n(&self, n: i32) -> Vec<f64> {
        assert!(n >= 0, "n must be non-negative, got {n}");
        let mut rng = crate::rng::single_rng();
        (0..n).map(|_| self.sample(&mut rng)).collect()
    }
}

// ── Family constructors (exposed to R as `dist_<name>`) ──────────────────────────
//
// The `dist_` prefix avoids masking base/stats functions (e.g. `base::gamma`,
// `stats::poisson`) when the package is attached.

/// Create a normal (Gaussian) distribution.
///
/// The second argument is the **variance** (σ²), not the standard deviation, to
/// match the way variance is usually quoted in statistical models. The standard
/// deviation passed to the underlying sampler is `sqrt(variance)`.
///
/// @param mean      Mean (μ) of the distribution.
/// @param variance  Variance (σ²); must be non-negative.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_normal(7, 4)   # mean 7, variance 4 (sd 2)
/// d$sample_one()
/// @export
#[extendr]
fn dist_normal(mean: f64, variance: f64) -> Distribution {
    assert!(
        variance >= 0.0 && variance.is_finite(),
        "variance must be finite and non-negative, got {variance}"
    );
    assert!(mean.is_finite(), "mean must be finite, got {mean}");
    let dist = Normal::new(mean, variance.sqrt())
        .unwrap_or_else(|e| panic!("invalid normal parameters (mean={mean}, variance={variance}): {e}"));
    Distribution {
        kind: DistKind::Normal(dist),
    }
}

/// Create a degenerate (constant) distribution that always returns `value`.
///
/// Use this as a fixed-duration drop-in wherever a `Distribution` is required —
/// e.g. a deterministic infectious period of exactly `value` ticks.
///
/// @param value The constant value returned by every draw.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_constant(10)   # always 10
/// d$sample_one()
/// @export
#[extendr]
fn dist_constant(value: f64) -> Distribution {
    assert!(value.is_finite(), "value must be finite, got {value}");
    Distribution {
        kind: DistKind::Constant(value),
    }
}

/// Create a continuous uniform distribution on the half-open interval [low, high).
///
/// Every value in `[low, high)` is equally likely. The mean is `(low + high) / 2`.
///
/// @param low   Inclusive lower bound.
/// @param high  Exclusive upper bound; must be strictly greater than `low`.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_uniform(3, 8)   # values in [3, 8), mean 5.5
/// d$sample_one()
/// @export
#[extendr]
fn dist_uniform(low: f64, high: f64) -> Distribution {
    assert!(
        low.is_finite() && high.is_finite(),
        "uniform bounds must be finite, got low={low}, high={high}"
    );
    assert!(
        high > low,
        "uniform requires high > low, got low={low}, high={high}"
    );
    Distribution {
        kind: DistKind::Uniform(Uniform::new(low, high)),
    }
}

/// Create a gamma distribution parameterized by shape and scale.
///
/// Uses the shape–scale (k, θ) parameterization: the mean is `shape * scale` and
/// the variance is `shape * scale^2`. Draws are strictly positive, which makes the
/// gamma a natural choice for right-skewed, always-positive durations.
///
/// @param shape  Shape parameter k; must be finite and positive.
/// @param scale  Scale parameter θ; must be finite and positive.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_gamma(2, 3)   # mean 6, variance 18
/// d$sample_one()
/// @export
#[extendr]
fn dist_gamma(shape: f64, scale: f64) -> Distribution {
    assert!(
        shape.is_finite() && shape > 0.0,
        "gamma shape must be finite and positive, got {shape}"
    );
    assert!(
        scale.is_finite() && scale > 0.0,
        "gamma scale must be finite and positive, got {scale}"
    );
    let dist = Gamma::new(shape, scale)
        .unwrap_or_else(|e| panic!("invalid gamma parameters (shape={shape}, scale={scale}): {e}"));
    Distribution {
        kind: DistKind::Gamma(dist),
    }
}

/// Create a Poisson distribution with rate (mean) `lambda`.
///
/// Draws are non-negative integer counts (returned as doubles) with mean and
/// variance both equal to `lambda`. Useful for count-valued durations.
///
/// @param lambda  Rate / mean λ; must be finite and positive.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_poisson(5)   # mean 5, integer-valued draws
/// d$sample_one()
/// @export
#[extendr]
fn dist_poisson(lambda: f64) -> Distribution {
    assert!(
        lambda.is_finite() && lambda > 0.0,
        "poisson lambda must be finite and positive, got {lambda}"
    );
    let dist = Poisson::new(lambda)
        .unwrap_or_else(|e| panic!("invalid poisson parameter (lambda={lambda}): {e}"));
    Distribution {
        kind: DistKind::Poisson(dist),
    }
}

/// Create a beta distribution on the open interval (0, 1).
///
/// Parameterized by the two positive shape parameters α and β. The mean is
/// `alpha / (alpha + beta)`.
///
/// @param alpha  First shape parameter (α); must be finite and positive.
/// @param beta   Second shape parameter (β); must be finite and positive.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_beta(2, 5)   # mean 2/7 ≈ 0.286, support (0, 1)
/// d$sample_one()
/// @export
#[extendr]
fn dist_beta(alpha: f64, beta: f64) -> Distribution {
    assert!(
        alpha.is_finite() && alpha > 0.0,
        "beta alpha must be finite and positive, got {alpha}"
    );
    assert!(
        beta.is_finite() && beta > 0.0,
        "beta beta must be finite and positive, got {beta}"
    );
    let dist = Beta::new(alpha, beta)
        .unwrap_or_else(|e| panic!("invalid beta parameters (alpha={alpha}, beta={beta}): {e}"));
    Distribution {
        kind: DistKind::Beta(dist),
    }
}

/// Create an exponential distribution with rate `rate`.
///
/// Draws are strictly positive with mean `1 / rate` and variance `1 / rate^2`.
///
/// @param rate  Rate parameter λ; must be finite and positive.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_exp(0.5)   # mean 2
/// d$sample_one()
/// @export
#[extendr]
fn dist_exp(rate: f64) -> Distribution {
    assert!(
        rate.is_finite() && rate > 0.0,
        "exponential rate must be finite and positive, got {rate}"
    );
    let dist = Exp::new(rate)
        .unwrap_or_else(|e| panic!("invalid exponential parameter (rate={rate}): {e}"));
    Distribution {
        kind: DistKind::Exp(dist),
    }
}

/// Create a logistic distribution with the given location and scale.
///
/// Symmetric about `location` (its mean and median); the variance is
/// `scale^2 * pi^2 / 3`.
///
/// @param location  Location parameter μ (the mean); must be finite.
/// @param scale     Scale parameter s; must be finite and positive.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_logistic(4, 2)   # mean 4
/// d$sample_one()
/// @export
#[extendr]
fn dist_logistic(location: f64, scale: f64) -> Distribution {
    assert!(location.is_finite(), "logistic location must be finite, got {location}");
    assert!(
        scale.is_finite() && scale > 0.0,
        "logistic scale must be finite and positive, got {scale}"
    );
    Distribution {
        kind: DistKind::Logistic { location, scale },
    }
}

/// Create a log-normal distribution.
///
/// A variable whose natural logarithm is `Normal(meanlog, sdlog)`. Draws are
/// strictly positive. The median is `exp(meanlog)` and the mean is
/// `exp(meanlog + sdlog^2 / 2)`. `meanlog` and `sdlog` are the log-space
/// parameters, matching R's `qlnorm(p, meanlog, sdlog)`.
///
/// @param meanlog  Mean of the underlying normal (in log space); must be finite.
/// @param sdlog    Standard deviation of the underlying normal; must be finite
///   and non-negative.
/// @return A `Distribution` object.
/// @examples
/// d <- dist_lognormal(0, 0.5)   # median 1
/// d$sample_one()
/// @export
#[extendr]
fn dist_lognormal(meanlog: f64, sdlog: f64) -> Distribution {
    assert!(meanlog.is_finite(), "lognormal meanlog must be finite, got {meanlog}");
    assert!(
        sdlog.is_finite() && sdlog >= 0.0,
        "lognormal sdlog must be finite and non-negative, got {sdlog}"
    );
    let dist = LogNormal::new(meanlog, sdlog)
        .unwrap_or_else(|e| panic!("invalid lognormal parameters (meanlog={meanlog}, sdlog={sdlog}): {e}"));
    Distribution {
        kind: DistKind::LogNormal(dist),
    }
}

extendr_module! {
    mod distributions;
    impl Distribution;
    fn dist_normal;
    fn dist_constant;
    fn dist_uniform;
    fn dist_gamma;
    fn dist_poisson;
    fn dist_beta;
    fn dist_exp;
    fn dist_logistic;
    fn dist_lognormal;
}
