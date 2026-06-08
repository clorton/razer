"""Tests for the worldwide agent-based SEIR model.

Failure here means the epidemic engine is wrong: bad FOI/R0, non-conservation of people,
broken timer progression, or spatial coupling that leaks or fails to spread.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.optimize import brentq

from wwsim.abm.model import SEIRParams, WorldSEIR
from wwsim.abm.networks import row_normalize


def _coupling(dense: np.ndarray):
    """Make a ``(WT, rowsum)`` coupling pair from a small dense weight matrix."""
    W = sp.csr_matrix(dense)
    return W.T.tocsr(), np.asarray(W.sum(axis=1)).ravel().astype(float)


def _nodes(pops):
    """A minimal node table (global_nodeid, population, iso3) for N nodes."""
    return pd.DataFrame({
        "global_nodeid": range(len(pops)),
        "population": pops,
        "iso3": ["AAA"] * len(pops),
        "adm2_name": [f"n{i}" for i in range(len(pops))],
    })


def test_row_normalize_rows_sum_to_fraction():
    """Given raw weights, when row-normalized, then each non-empty row sums to the fraction."""
    W = sp.csr_matrix(np.array([[0.0, 2.0, 6.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    Wn = row_normalize(W, 0.1)
    rowsum = np.asarray(Wn.sum(axis=1)).ravel()
    assert rowsum[0] == pytest.approx(0.1)
    assert rowsum[1] == 0.0          # empty row stays empty
    assert rowsum[2] == pytest.approx(0.1)


def test_population_is_conserved_every_tick():
    """Given a closed model, when run, then S+E+I+R equals the agent count at every tick.

    Agents only change disease state (never node), so the population must be invariant.
    """
    nodes = _nodes([3000, 2000, 1000])
    intra = _coupling(np.array([[0, .05, .05], [.05, 0, .05], [.05, .05, 0]]))
    p = SEIRParams(nticks=40, beta=0.4, incubation_period=3, infectious_period=5,
                   subsample=1, seed_count=20, seed_nodeid=0)
    m = WorldSEIR(nodes, p, intra)
    total = int(nodes["population"].sum())
    m.run()
    per_tick = sum(getattr(m.nodes, s) for s in ("S", "E", "I", "R")).sum(axis=1)
    assert np.all(per_tick == total)


def test_timer_progression_seir_timing():
    """Given seeded infectious and beta=0, when run, then they recover after infectious_period.

    With no transmission, the only dynamics are the timer: I -> R after `infectious_period`
    days, and nobody is ever exposed.
    """
    nodes = _nodes([10_000])
    intra = _coupling(np.zeros((1, 1)))
    D = 5
    p = SEIRParams(nticks=20, beta=0.0, incubation_period=3, infectious_period=D,
                   subsample=1, seed_count=200, seed_nodeid=0)
    m = WorldSEIR(nodes, p, intra)
    m.run()
    tot = m.totals()
    assert tot["E"].max() == 0                       # beta=0 => no exposures
    assert tot["I"].iloc[0] == 200                   # seeded infectious present at start
    assert tot["I"].iloc[D + 1] == 0                 # all recovered after the infectious period
    assert tot["R"].iloc[-1] == 200                  # exactly the seeds recovered


def test_intra_coupling_spreads_to_neighbor():
    """Given two coupled nodes, when one is seeded, then infection reaches the other.

    The sparse FOI coupling must export infection pressure along network edges.
    """
    nodes = _nodes([20_000, 20_000])
    intra = _coupling(np.array([[0.0, 0.1], [0.1, 0.0]]))
    p = SEIRParams(nticks=120, beta=0.5, incubation_period=3, infectious_period=6,
                   subsample=1, seed_count=50, seed_nodeid=0)
    m = WorldSEIR(nodes, p, intra)
    m.run()
    assert m.nodes.R[-1, 1] > 0  # node 1 (no seeds) became infected via coupling


def test_no_coupling_isolates_nodes():
    """Given two uncoupled nodes, when one is seeded, then the other never gets infected.

    Guards against an FOI bug that mixes nodes regardless of the network.
    """
    nodes = _nodes([20_000, 20_000])
    intra = _coupling(np.zeros((2, 2)))
    p = SEIRParams(nticks=120, beta=0.5, incubation_period=3, infectious_period=6,
                   subsample=1, seed_count=50, seed_nodeid=0)
    m = WorldSEIR(nodes, p, intra)
    m.run()
    assert m.nodes.R[-1, 1] == 0  # node 1 stays untouched


def test_forcing_unity_equals_unforced():
    """Given a forcing of all 1s, when run, then it matches the unforced run exactly.

    A multiplier of 1 must be a no-op; any difference means the forcing path corrupts beta.
    """
    nodes = _nodes([20_000])
    intra = _coupling(np.zeros((1, 1)))
    p = SEIRParams(nticks=120, beta=0.5, incubation_period=3, infectious_period=6,
                   subsample=1, seed_count=50, seed_nodeid=0)
    base = WorldSEIR(nodes, p, intra)
    base.run()
    forced = WorldSEIR(nodes, p, intra, forcing=(np.ones((120, 1), np.float32), ["AAA"]))
    forced.run()
    assert forced.nodes.R[-1, 0] == base.nodes.R[-1, 0]


def test_forcing_zero_suppresses_transmission():
    """Given a forcing of 0, when run, then nobody new is infected (only seeds recover).

    m=0 zeroes beta, so the country's FOI is 0 -- a hard check that the multiplier gates FOI.
    """
    nodes = _nodes([20_000])
    intra = _coupling(np.zeros((1, 1)))
    p = SEIRParams(nticks=120, beta=0.5, incubation_period=3, infectious_period=6,
                   subsample=1, seed_count=50, seed_nodeid=0)
    m = WorldSEIR(nodes, p, intra, forcing=(np.zeros((120, 1), np.float32), ["AAA"]))
    m.run()
    assert m.totals()["E"].max() == 0          # no exposures ever
    assert m.nodes.R[-1, 0] == 50              # exactly the seeds recovered


def test_forcing_above_one_increases_attack():
    """Given a forcing > 1, when run, then the attack fraction exceeds the unforced run."""
    nodes = _nodes([20_000])
    intra = _coupling(np.zeros((1, 1)))
    p = SEIRParams(nticks=160, beta=0.25, incubation_period=3, infectious_period=6,
                   subsample=1, seed_count=50, seed_nodeid=0)
    base = WorldSEIR(nodes, p, intra)
    base.run()
    hotter = WorldSEIR(nodes, p, intra, forcing=(np.full((160, 1), 1.6, np.float32), ["AAA"]))
    hotter.run()
    assert hotter.nodes.R[-1, 0] > base.nodes.R[-1, 0]


@pytest.mark.slow
def test_attack_fraction_matches_kermack_mckendrick():
    """Given a single well-mixed node, when run, then attack fraction ≈ KM final size.

    The core epidemiological correctness check: R0 = beta * infectious_period and the final
    size solves 1 - z = exp(-R0 z). A mismatch means the FOI or progression is wrong.
    """
    N = 120_000
    nodes = _nodes([N])
    intra = _coupling(np.zeros((1, 1)))
    beta, D = 0.3, 6
    p = SEIRParams(nticks=300, beta=beta, incubation_period=4, infectious_period=D,
                   subsample=1, seed_count=50, seed_nodeid=0)
    m = WorldSEIR(nodes, p, intra)
    m.run()
    attack = m.nodes.R[-1, 0] / N
    R0 = beta * D
    z = brentq(lambda z: 1 - z - np.exp(-R0 * z), 1e-9, 1 - 1e-9)
    assert attack == pytest.approx(z, abs=0.02)
