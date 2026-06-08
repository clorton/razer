"""Build the one global cross-border air-travel network over admin-2 nodes.

Pipeline (the goal, restated):

1. Start from weighted directed **airport** edges (:func:`wwsim.flights.build_airport_edges`).
2. Keep only the **top-N airports by passenger volume** (both endpoints must be top-N).
3. Map each airport to its **admin-2 node** (:func:`wwsim.airports_assign.assign_airports_to_nodes`).
4. Keep only **cross-border** edges -- the two airports lie in different countries.
5. **Aggregate** all airports in the same admin-2 unit: sum edge weights to get directed
   admin-2-to-admin-2 edges. So e.g. all of London's airports collapse into one London
   admin-2 node, and an edge London->New-York carries the combined traffic of every
   airport-pair between those two admin-2 units.

The result is one global network whose nodes are admin-2 units and whose edges are
cross-border air-travel volumes -- the inter-nation coupling for the worldwide model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import Config
from .flights import airport_volume
from .logging import logger


def build_air_network(
    airport_edges: pd.DataFrame,
    airport_nodes: pd.DataFrame,
    cfg: Config,
) -> pd.DataFrame:
    """Aggregate airport edges into a cross-border admin-2-to-admin-2 air network.

    Args:
        airport_edges: Directed weighted airport edges from
            :func:`wwsim.flights.build_airport_edges` (``src_iata, dst_iata, weight,
            n_carriers``).
        airport_nodes: Airport-to-node assignment from
            :func:`wwsim.airports_assign.assign_airports_to_nodes` (``iata,
            node_global_nodeid, node_iso3``).
        cfg: Project configuration (uses ``cfg.air``: ``top_n_airports``,
            ``intra_country_air``).

    Returns:
        DataFrame ``[src_global_nodeid, dst_global_nodeid, src_iso3, dst_iso3, weight,
        n_airport_pairs, n_carriers]`` -- directed admin-2 edges, cross-border only (unless
        ``cfg.air.intra_country_air``).
    """
    ap = cfg.air

    # --- top-N airports by passenger-proxy volume ---
    edges = airport_edges
    if ap.top_n_airports is not None:
        vol = airport_volume(edges)
        keep = set(vol.head(ap.top_n_airports).index)
        before = len(edges)
        edges = edges[edges["src_iata"].isin(keep) & edges["dst_iata"].isin(keep)]
        logger.info(
            "air_network: top-%d airports -> %d/%d airport edges retained",
            ap.top_n_airports, len(edges), before,
        )

    # --- map airports to nodes ---
    amap = airport_nodes.set_index("iata")[["node_global_nodeid", "node_iso3"]]
    edges = edges.join(amap.rename(columns={"node_global_nodeid": "src_node", "node_iso3": "src_node_iso3"}), on="src_iata")
    edges = edges.join(amap.rename(columns={"node_global_nodeid": "dst_node", "node_iso3": "dst_node_iso3"}), on="dst_iata")
    edges = edges.dropna(subset=["src_node", "dst_node", "src_node_iso3", "dst_node_iso3"])
    edges["src_node"] = edges["src_node"].astype("int64")
    edges["dst_node"] = edges["dst_node"].astype("int64")

    # --- cross-border filter (the global network carries only between-country air travel) ---
    if not ap.intra_country_air:
        before = len(edges)
        edges = edges[edges["src_node_iso3"] != edges["dst_node_iso3"]]
        logger.info("air_network: cross-border filter -> %d/%d edges", len(edges), before)

    # Drop self-loops that can arise when two airports share one admin-2 node.
    edges = edges[edges["src_node"] != edges["dst_node"]]

    # --- aggregate airports within the same admin-2 node ---
    agg = (
        edges.groupby(["src_node", "dst_node"], as_index=False)
        .agg(
            src_iso3=("src_node_iso3", "first"),
            dst_iso3=("dst_node_iso3", "first"),
            weight=("weight", "sum"),
            n_airport_pairs=("weight", "size"),
            n_carriers=("n_carriers", "sum"),
        )
        .rename(columns={"src_node": "src_global_nodeid", "dst_node": "dst_global_nodeid"})
    )
    logger.info(
        "air_network: %d directed admin-2 edges across %d source / %d dest nodes",
        len(agg), agg["src_global_nodeid"].nunique(), agg["dst_global_nodeid"].nunique(),
    )
    return agg


def air_network_to_sparse(edges: pd.DataFrame, n_nodes: int) -> sp.csr_matrix:
    """Convert admin-2 air edges to a global sparse adjacency matrix.

    Args:
        edges: Output of :func:`build_air_network`.
        n_nodes: Total number of global nodes (matrix dimension).

    Returns:
        An ``(n_nodes, n_nodes)`` CSR matrix of air-travel weights.
    """
    if len(edges) == 0:
        return sp.csr_matrix((n_nodes, n_nodes))
    return sp.coo_matrix(
        (
            edges["weight"].to_numpy(),
            (edges["src_global_nodeid"].to_numpy(), edges["dst_global_nodeid"].to_numpy()),
        ),
        shape=(n_nodes, n_nodes),
    ).tocsr()


def save_air_network(
    edges: pd.DataFrame, cfg: Config, n_nodes: int, stem: str = "global_air_network"
) -> None:
    """Persist the air network as a Parquet edge list and a sparse ``.npz`` matrix.

    Args:
        edges: Output of :func:`build_air_network`.
        cfg: Project configuration.
        n_nodes: Total number of global nodes.
        stem: Output filename stem (default ``"global_air_network"``); top-N variants use
            e.g. ``"global_air_network_top250"`` so several can coexist.
    """
    cfg.networks_dir.mkdir(parents=True, exist_ok=True)
    edges.to_parquet(cfg.networks_dir / f"{stem}.parquet", index=False)
    sp.save_npz(cfg.networks_dir / f"{stem}.npz", air_network_to_sparse(edges, n_nodes))
    logger.info("air_network: saved %s.{parquet,npz}", stem)
