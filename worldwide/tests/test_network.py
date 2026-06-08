"""Tests for the pluggable multi-modal combiner.

Failure here means modes (gravity, air, rail) are not summed correctly into the global
matrix, so the simulation's coupling is wrong.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from wwsim.network import ModeNetwork, combine_modes, gravity_edges_from_matrices


def test_combine_modes_sums_overlapping_edges_with_scale():
    """Given two modes sharing an edge, when combined, then scaled weights add at that cell.

    This is the core of multi-modality: air + rail on the same node pair must accumulate.
    """
    gravity = ModeNetwork(
        "gravity",
        pd.DataFrame({"src_global_nodeid": [0], "dst_global_nodeid": [1], "weight": [2.0]}),
        scale=1.0,
    )
    air = ModeNetwork(
        "air",
        pd.DataFrame({"src_global_nodeid": [0], "dst_global_nodeid": [1], "weight": [10.0]}),
        scale=0.5,
    )
    m = combine_modes([gravity, air], n_nodes=2)
    assert m[0, 1] == 2.0 + 0.5 * 10.0  # 2 + 5 = 7
    assert m.shape == (2, 2)


def test_combine_modes_empty_is_zero_matrix():
    """Given no edges, when combined, then an all-zero matrix of the right size results."""
    m = combine_modes([], n_nodes=3)
    assert m.shape == (3, 3)
    assert m.nnz == 0


def test_gravity_edges_from_matrices_flattens_all_countries():
    """Given per-country matrices, when flattened, then all non-zero edges are concatenated."""
    matrices = {
        "AAA": (np.array([0, 1]), np.array([[0.0, 3.0], [4.0, 0.0]])),
        "BBB": (np.array([2, 3]), np.array([[0.0, 0.0], [5.0, 0.0]])),
    }
    edges = gravity_edges_from_matrices(matrices)
    assert len(edges) == 3  # 2 from AAA + 1 from BBB
    assert set(edges["src_global_nodeid"]) == {0, 1, 3}
