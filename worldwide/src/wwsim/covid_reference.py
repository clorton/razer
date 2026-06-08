"""Real COVID-19 reference data (Our World in Data) for comparison with the sim.

Note on dates: COVID-19 (SARS-CoV-2) emerged in late 2019; the pandemic ran ~2020-2024.
There is no COVID data for 2000-2004. This module fetches Our World in Data's per-country
daily cases/deaths (the standard open aggregate, compiled from WHO/JHU) and plots the global
and per-country curves so they can be compared, qualitatively, to a :mod:`wwsim.abm` run.

What is comparable, and what is not: our SEIR run is a single, intervention-free wave of
*true infections* on 2015 populations seeded in one place. The real series are *reported
cases* (shaped by testing capacity, NPIs/lockdowns, and successive variants -- Alpha, Delta,
Omicron), so expect the real world to show several waves and far more structure. The honest
comparison is the *shape of a single wave* and the relative ordering of country arrivals.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pycountry  # noqa: E402
import requests  # noqa: E402

from .config import Config  # noqa: E402
from .logging import logger  # noqa: E402

# OWID location names that pycountry's lookup does not resolve to the right alpha-3.
_OWID_ISO_OVERRIDES = {
    "Democratic Republic of Congo": "COD", "Congo": "COG", "Cote d'Ivoire": "CIV",
    "South Korea": "KOR", "North Korea": "PRK", "Laos": "LAO", "Brunei": "BRN",
    "Cape Verde": "CPV", "Cabo Verde": "CPV", "East Timor": "TLS", "Timor": "TLS",
    "Micronesia (country)": "FSM", "Sao Tome and Principe": "STP", "Syria": "SYR",
    "Russia": "RUS", "Iran": "IRN", "Vietnam": "VNM", "Moldova": "MDA", "Bolivia": "BOL",
    "Venezuela": "VEN", "Tanzania": "TZA", "Turkey": "TUR", "Turkiye": "TUR",
    "United States": "USA", "United Kingdom": "GBR", "Czechia": "CZE", "Eswatini": "SWZ",
    "Gambia": "GMB", "Kyrgyzstan": "KGZ",
}


def owid_to_iso3(name: str) -> str | None:
    """Map an OWID location name to an ISO 3166-1 alpha-3 code (or ``None``).

    Aggregates (``World``, continents, income groups) and non-UN entities return ``None``.

    Args:
        name: OWID ``location`` value.

    Returns:
        Alpha-3 code or ``None`` if it is not a single mappable country.
    """
    if name in _OWID_ISO_OVERRIDES:
        return _OWID_ISO_OVERRIDES[name]
    try:
        return pycountry.countries.lookup(name).alpha_3
    except LookupError:
        return None


def country_value_matrix(
    df: pd.DataFrame,
    iso3_order: list[str],
    populations: dict[str, float],
    n_days: int,
    field: str = "new_cases",
    per_capita_per: float = 100_000,
    smooth: int = 7,
    start_date: pd.Timestamp | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a dense ``(n_days, n_countries)`` per-capita, smoothed case/death matrix.

    Args:
        df: Output of :func:`load_owid`.
        iso3_order: Column order (the countries to include, e.g. the 193 UN members).
        populations: ISO3 -> population, for the per-capita rate.
        n_days: Number of days from ``start_date``.
        field: ``new_cases`` or ``new_deaths``.
        per_capita_per: Rate denominator (per 100,000 by default).
        smooth: Rolling-mean window (days).
        start_date: First day; defaults to the dataset's first date.

    Returns:
        Tuple ``(dates[n_days] datetime64, values[n_days, n_countries] float32)``. Countries
        with no OWID data are all-``NaN`` columns (rendered as "no data").
    """
    d = df.copy()
    d["iso3"] = d["location"].map(owid_to_iso3)
    d = d[d["iso3"].isin(set(iso3_order))]

    start = start_date if start_date is not None else d["date"].min()
    dates = pd.date_range(start, periods=n_days, freq="D")

    wide = (
        d.pivot_table(index="date", columns="iso3", values=field, aggfunc="sum")
        .reindex(dates)
        .reindex(columns=iso3_order)
    )
    wide = wide.clip(lower=0).rolling(smooth, min_periods=1).mean()

    pop = np.array([populations.get(c, np.nan) for c in iso3_order], dtype=np.float64)
    values = wide.to_numpy(dtype=np.float64) / pop[None, :] * per_capita_per
    n_mapped = int(np.isfinite(values).any(axis=0).sum())
    logger.info("covid: country matrix %dx%d, %d/%d countries have data, field=%s",
                n_days, len(iso3_order), n_mapped, len(iso3_order), field)
    return dates.to_numpy(), values.astype(np.float32)

OWID_URL = (
    "https://raw.githubusercontent.com/owid/covid-19-data/master/"
    "public/data/cases_deaths/full_data.csv"
)
# The mega-file carries the per-country effective reproduction number (`reproduction_rate`,
# an Arroyo-Marioli et al. Kalman-filter estimate from daily cases), keyed by ISO3.
OWID_FULL_URL = (
    "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"
)


def download_owid(cfg: Config, force: bool = False) -> Path:
    """Download the OWID cases/deaths CSV into ``data/covid/`` (cached).

    Args:
        cfg: Project configuration.
        force: Re-download even if present.

    Returns:
        Path to the local CSV.

    Raises:
        requests.HTTPError: If the download fails.
    """
    dest = cfg.data_dir / "covid" / "owid_full_data.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        logger.info("covid: using cached %s", dest.name)
        return dest
    logger.info("covid: downloading OWID cases/deaths...")
    resp = requests.get(OWID_URL, timeout=180)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return dest


def download_owid_full(cfg: Config, force: bool = False) -> Path:
    """Download the OWID mega-file (with ``reproduction_rate``) into ``data/covid/`` (cached).

    Args:
        cfg: Project configuration.
        force: Re-download even if present.

    Returns:
        Path to the local CSV (~94 MB).

    Raises:
        requests.HTTPError: If the download fails.
    """
    dest = cfg.data_dir / "covid" / "owid_covid_data.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        logger.info("covid: using cached %s", dest.name)
        return dest
    logger.info("covid: downloading OWID mega-file (~94 MB) for reproduction_rate...")
    with requests.get(OWID_FULL_URL, timeout=300, stream=True) as resp:
        resp.raise_for_status()
        with dest.open("wb") as f:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    return dest


def country_rt_matrix(cfg: Config, iso3_order: list[str], dates: pd.DatetimeIndex) -> np.ndarray:
    """Per-country effective reproduction number aligned to a calendar date axis.

    Args:
        cfg: Project configuration.
        iso3_order: Column order (ISO3 codes).
        dates: Calendar dates (rows) to align Rt to.

    Returns:
        ``(len(dates), n_countries)`` float64 array of Rt, ``NaN`` where OWID has no estimate.
        Small internal gaps (<= 7 days) are forward-filled within a country.
    """
    path = download_owid_full(cfg)
    df = pd.read_csv(path, usecols=["iso_code", "date", "reproduction_rate"], parse_dates=["date"])
    df = df[df["iso_code"].isin(set(iso3_order))]
    wide = (
        df.pivot_table(index="date", columns="iso_code", values="reproduction_rate", aggfunc="mean")
        .reindex(dates)
        .ffill(limit=7)            # bridge short reporting gaps; leave long gaps as NaN
        .reindex(columns=iso3_order)
    )
    return wide.to_numpy(dtype=np.float64)


def rt_forcing(
    cfg: Config,
    iso3_order: list[str],
    anchor_date: pd.Timestamp,
    nticks: int,
    baseline: float,
    clip_max: float = 5.0,
) -> np.ndarray:
    """Build the per-country, per-tick transmission multiplier ``m_c(t) = Rt_c(t) / baseline``.

    Sim tick ``t`` maps to calendar date ``anchor_date + t days``. The multiplier scales the
    model's ``beta`` per country, so with ``baseline = R0`` the local instantaneous
    reproduction number tracks the empirical Rt (the model still applies susceptible
    depletion mechanistically on top). Missing Rt → ``m = 1`` (no forcing).

    Args:
        cfg: Project configuration.
        iso3_order: Country column order (must match the model's mapping).
        anchor_date: Calendar date for sim day 0.
        nticks: Number of sim ticks (days).
        baseline: Divisor (use the model's ``R0 = beta * infectious_period`` for m=1 at Rt=R0).
        clip_max: Upper clamp on the multiplier (guards against noisy Rt spikes).

    Returns:
        ``(nticks, n_countries)`` float32 multiplier array.
    """
    dates = pd.date_range(anchor_date, periods=nticks, freq="D")
    rt = country_rt_matrix(cfg, iso3_order, dates)            # (nticks, C), NaN where missing
    m = rt / baseline
    m = np.where(np.isfinite(m), m, 1.0)                      # no data -> no forcing
    m = np.clip(m, 0.0, clip_max)
    coverage = float(np.isfinite(rt).mean())
    logger.info("covid: Rt forcing %dx%d, baseline=%.3g, real-Rt coverage=%.0f%% (rest m=1), "
                "mean m=%.3f", nticks, len(iso3_order), baseline, 100 * coverage, float(m.mean()))
    return m.astype(np.float32)


def load_owid(cfg: Config) -> pd.DataFrame:
    """Load the OWID CSV, parsing dates.

    Args:
        cfg: Project configuration.

    Returns:
        Long-format DataFrame with ``date, location, new_cases, new_deaths, total_cases,
        total_deaths`` (and weekly/biweekly columns).
    """
    return pd.read_csv(cfg.data_dir / "covid" / "owid_full_data.csv", parse_dates=["date"])


def location_series(df: pd.DataFrame, location: str, smooth: int = 7) -> pd.DataFrame:
    """Per-location daily series with a rolling-mean smoothing of new cases/deaths.

    Args:
        df: Output of :func:`load_owid`.
        location: OWID location name (e.g. ``"World"``, ``"United States"``).
        smooth: Rolling-mean window in days for the smoothed columns.

    Returns:
        DataFrame indexed by date with ``new_cases, new_deaths, total_cases, total_deaths``
        and smoothed ``new_cases_smooth, new_deaths_smooth`` (negatives from reporting
        corrections clipped to 0).
    """
    s = df[df["location"] == location].sort_values("date").set_index("date")
    out = s[["new_cases", "new_deaths", "total_cases", "total_deaths"]].copy()
    out["new_cases"] = out["new_cases"].clip(lower=0)
    out["new_deaths"] = out["new_deaths"].clip(lower=0)
    out["new_cases_smooth"] = out["new_cases"].rolling(smooth, min_periods=1).mean()
    out["new_deaths_smooth"] = out["new_deaths"].rolling(smooth, min_periods=1).mean()
    return out


def plot_global(df: pd.DataFrame, out_path: Path) -> None:
    """Global daily new cases (7-day) + cumulative, with daily deaths on a twin axis.

    Args:
        df: Output of :func:`load_owid`.
        out_path: PNG path.
    """
    w = location_series(df, "World")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9), sharex=True)

    ax1.plot(w.index, w["new_cases_smooth"], color="crimson", lw=1.5, label="daily new cases (7-day avg)")
    ax1.set_ylabel("daily new cases", color="crimson")
    axc = ax1.twinx()
    axc.plot(w.index, w["total_cases"] / 1e9, color="navy", lw=1.5, label="cumulative cases")
    axc.set_ylabel("cumulative cases (billions)", color="navy")
    ax1.set_title(f"Global reported COVID-19 cases ({w.index.min().date()} to {w.index.max().date()}, OWID)")

    ax2.plot(w.index, w["new_deaths_smooth"], color="black", lw=1.5)
    ax2.set_ylabel("daily new deaths (7-day avg)")
    ax2.set_xlabel("date")
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    plt.close(fig)
    logger.info("covid: wrote %s", out_path)


def plot_countries(df: pd.DataFrame, locations: list[str], out_path: Path) -> None:
    """Per-country daily new cases (7-day avg) -- shows wave timing/ordering.

    Args:
        df: Output of :func:`load_owid`.
        locations: OWID location names.
        out_path: PNG path.
    """
    fig, ax = plt.subplots(figsize=(13, 7))
    for loc in locations:
        s = location_series(df, loc)
        if len(s):
            ax.plot(s.index, s["new_cases_smooth"], lw=1.5, label=loc)
    ax.set_xlabel("date")
    ax.set_ylabel("daily new cases (7-day avg)")
    ax.set_title("Reported COVID-19 cases by country (OWID)")
    ax.legend(ncol=2, fontsize=8)
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    plt.close(fig)
    logger.info("covid: wrote %s", out_path)


def plot_first_wave_vs_sim(df: pd.DataFrame, sim_totals: pd.DataFrame, out_path: Path,
                           wave_days: int = 200) -> None:
    """Overlay the shape of the real first global wave with the sim's incidence (normalized).

    Both curves are normalized to their own peak and aligned at day 0 (real day 0 = first day
    global cases exceed 1% of the first-wave peak), so only the *shape* is compared.

    Args:
        df: Output of :func:`load_owid`.
        sim_totals: A :meth:`wwsim.abm.model.WorldSEIR.totals` frame (uses ``incidence``).
        out_path: PNG path.
        wave_days: Number of days of the real series to treat as the first wave window.
    """
    w = location_series(df, "World")["new_cases_smooth"]
    first = w.iloc[:wave_days]
    start = (first > 0.01 * first.max()).idxmax()
    real = first.loc[start:].reset_index(drop=True)
    real_n = real / real.max()

    sim = sim_totals["incidence"].reset_index(drop=True)
    sim_n = sim / sim.max()

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(real_n.index, real_n.values, color="crimson", lw=2, label="real world first wave (2020)")
    ax.plot(sim_n.index, sim_n.values, color="steelblue", lw=2, label="wwsim global incidence")
    ax.set_xlabel("days since wave onset")
    ax.set_ylabel("new infections / cases (normalized to peak)")
    ax.set_title("Single-wave shape: real first wave vs simulation (normalized)")
    ax.legend()
    fig.text(0.5, -0.02,
             "Caveat: real = reported cases shaped by testing + NPIs; sim = intervention-free "
             "true infections on 2015 population. Shapes, not magnitudes, are comparable.",
             ha="center", fontsize=8, style="italic")
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    plt.close(fig)
    logger.info("covid: wrote %s", out_path)
