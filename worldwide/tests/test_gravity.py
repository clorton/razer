"""Tests for the gravity model.

Failure here means per-country migration matrices are wrong, distorting within-country
spatial spread in the simulation.
"""

from __future__ import annotations

import numpy as np

from wwsim.config import GravityParams
from wwsim.gravity import (
    edge_list_from_matrix,
    gravity_matrix,
    haversine_matrix,
    row_normalize,
)


def test_haversine_one_degree_latitude_is_about_111km():
    """Given two points 1 deg apart in latitude, when measured, then ~111 km results.

    One degree of latitude is ~111.2 km everywhere; a wrong constant breaks all distances.
    """
    d = haversine_matrix(np.array([0.0, 0.0]), np.array([0.0, 1.0]))
    assert abs(d[0, 1] - 111.19) < 1.0
    assert d[0, 0] == 0.0


def test_gravity_zero_diagonal_and_symmetry():
    """Given equal source/dest exponents, when built, then the matrix is symmetric with 0 diagonal."""
    params = GravityParams(k=1.0, a=1.0, b=1.0, c=2.0)
    pop = np.array([100.0, 200.0, 300.0])
    lon = np.array([0.0, 1.0, 2.0])
    lat = np.array([0.0, 0.0, 0.0])
    m = gravity_matrix(pop, lon, lat, params)
    assert np.allclose(np.diag(m), 0.0)
    assert np.allclose(m, m.T)  # a == b => symmetric


def test_gravity_decreases_with_distance():
    """Given a near and a far destination of equal pop, when built, then nearer is stronger.

    The 1/D^c decay is the core of the model; a sign/way error would invert spread.
    """
    params = GravityParams(k=1.0, a=1.0, b=1.0, c=2.0, min_distance_km=1.0)
    pop = np.array([100.0, 100.0, 100.0])
    lon = np.array([0.0, 1.0, 5.0])  # node1 near node0, node2 far
    lat = np.array([0.0, 0.0, 0.0])
    m = gravity_matrix(pop, lon, lat, params)
    assert m[0, 1] > m[0, 2]


def test_row_normalize_rows_sum_to_one():
    """Given a matrix, when row-normalized, then non-empty rows sum to 1."""
    m = np.array([[0.0, 1.0, 3.0], [0.0, 0.0, 0.0], [2.0, 2.0, 0.0]])
    n = row_normalize(m)
    assert abs(n[0].sum() - 1.0) < 1e-12
    assert n[1].sum() == 0.0  # all-zero row stays zero
    assert abs(n[2].sum() - 1.0) < 1e-12


def test_edge_list_only_nonzero():
    """Given a matrix, when flattened to edges, then only non-zero entries appear."""
    ids = np.array([10, 20])
    m = np.array([[0.0, 5.0], [0.0, 0.0]])
    edges = edge_list_from_matrix(ids, m)
    assert len(edges) == 1
    assert edges.iloc[0]["src_global_nodeid"] == 10
    assert edges.iloc[0]["dst_global_nodeid"] == 20
    assert edges.iloc[0]["weight"] == 5.0
