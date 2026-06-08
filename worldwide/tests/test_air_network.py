"""Tests for the global cross-border air network.

Failure here means the inter-country coupling is wrong: intra-country edges leak in, the
top-N filter misfires, or airports in one admin-2 are not aggregated.
"""

from __future__ import annotations

import pandas as pd

from wwsim.air_network import air_network_to_sparse, build_air_network
from wwsim.config import Config


def _edges():
    """Airport edges: two AAA airports -> one BBB airport, plus an intra-AAA edge."""
    return pd.DataFrame(
        {
            "src_iata": ["AA1", "AA2", "AA1"],
            "dst_iata": ["BB1", "BB1", "AA2"],
            "weight": [100.0, 50.0, 999.0],
            "n_carriers": [2, 1, 5],
            "src_iso3": ["AAA", "AAA", "AAA"],
            "dst_iso3": ["BBB", "BBB", "AAA"],
        }
    )


def _assignment():
    """AA1 and AA2 both sit in AAA node 0; BB1 in BBB node 2."""
    return pd.DataFrame(
        {
            "iata": ["AA1", "AA2", "BB1"],
            "node_global_nodeid": [0, 0, 2],
            "node_iso3": ["AAA", "AAA", "BBB"],
        }
    )


def test_cross_border_only_and_aggregated():
    """Given mixed edges, when built, then only cross-border survives, aggregated per admin-2.

    The intra-AAA edge (AA1->AA2) must be dropped; the two AAA->BBB airport edges must
    collapse into one node-0 -> node-2 edge with summed weight.
    """
    cfg = Config()
    cfg.air.top_n_airports = None
    air = build_air_network(_edges(), _assignment(), cfg)
    assert len(air) == 1
    row = air.iloc[0]
    assert row["src_global_nodeid"] == 0 and row["dst_global_nodeid"] == 2
    assert row["src_iso3"] == "AAA" and row["dst_iso3"] == "BBB"
    assert row["weight"] == 150.0  # 100 + 50 aggregated
    assert row["n_airport_pairs"] == 2


def test_top_n_airport_filter_drops_small_airports():
    """Given a top-N of 2, when built, then routes touching the 3rd-busiest airport drop.

    Volumes: AA1=100+999=1099, AA2=50+999=1049, BB1=150. top-2 = {AA1, AA2}, so the only
    surviving edges need both endpoints in {AA1, AA2}; the cross-border AAA->BBB edges use
    BB1 and are removed -> empty cross-border network.
    """
    cfg = Config()
    cfg.air.top_n_airports = 2
    air = build_air_network(_edges(), _assignment(), cfg)
    assert len(air) == 0


def test_to_sparse_places_weight_at_node_indices():
    """Given air edges, when converted to sparse, then weight lands at [src, dst]."""
    cfg = Config()
    cfg.air.top_n_airports = None
    air = build_air_network(_edges(), _assignment(), cfg)
    mat = air_network_to_sparse(air, n_nodes=4)
    assert mat[0, 2] == 150.0
    assert mat.shape == (4, 4)
    assert mat.nnz == 1
