"""Sparse migration networks for the ABM force-of-infection coupling.

laser-generic couples nodes with a **dense** ``ft[:, None] * network`` operation. At 45,406
nodes that matrix is ~8 GB and ~2 billion multiplies per tick -- infeasible. The dynamics it
encodes, though, are sparse and cheap: a node exports a fraction of its force of infection
(FOI) to the nodes it is connected to. This module assembles the two coupling layers as
**sparse** matrices and pre-transposes them so the per-tick FOI step is two sparse
matrix-vector products (see :class:`wwsim.abm.components.Transmission`).

- **Intra-country** coupling: each country's gravity matrix, placed block-diagonally over
  the global node index. Within-country spread only.
- **Air** coupling (optional): the global cross-border air network for the chosen top-N
  airports. Only admin-2 nodes that contain a selected airport have nonzero rows/cols, so
  the coupling touches exactly the gateway nodes -- and only across borders.

Each layer is **row-normalized to a fixed out-fraction**: every node with out-edges exports
that fraction of its FOI, split among destinations by the raw (gravity / seat) weights. The
fraction is the knob for how strongly space couples the epidemic.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import Config
from ..gravity import build_all_country_matrices, country_matrices_to_global_sparse
from ..logging import logger


def row_normalize(W: sp.csr_matrix, out_fraction: float) -> sp.csr_matrix:
    """Scale each row so it sums to ``out_fraction`` (rows that sum to 0 stay 0).

    Args:
        W: Raw weighted adjacency (CSR), nonnegative.
        out_fraction: Target row sum -- the fraction of a node's FOI it exports.

    Returns:
        Row-normalized CSR matrix.
    """
    rowsum = np.asarray(W.sum(axis=1)).ravel()
    inv = np.zeros_like(rowsum, dtype=np.float64)
    nz = rowsum > 0
    inv[nz] = out_fraction / rowsum[nz]
    return (sp.diags(inv) @ W).tocsr()


def _transpose_and_rowsum(Wn: sp.csr_matrix) -> tuple[sp.csr_matrix, np.ndarray]:
    """Return (transpose as CSR, row sums) -- the two things the FOI step needs.

    The transpose gives **incoming** pressure (``Wᵀ @ ft``); the row sums give **outgoing**
    pressure (``ft * rowsum``). Computing both once avoids per-tick work.
    """
    rowsum = np.asarray(Wn.sum(axis=1)).ravel().astype(np.float64)
    return Wn.T.tocsr(), rowsum


def build_intra_coupling(
    global_nodes: pd.DataFrame, cfg: Config, out_fraction: float
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Assemble the block-diagonal per-country gravity coupling.

    Args:
        global_nodes: The global node table (defines node count and per-country grouping).
        cfg: Project configuration (gravity params).
        out_fraction: FOI export fraction per node.

    Returns:
        Tuple ``(WT, rowsum)`` where ``WT`` is the transposed row-normalized network (CSR)
        and ``rowsum`` is the per-node exported fraction.
    """
    matrices = build_all_country_matrices(global_nodes, cfg, save=False)
    W = country_matrices_to_global_sparse(matrices, len(global_nodes))
    Wn = row_normalize(W, out_fraction)
    logger.info("abm.networks: intra coupling %d nodes, %d edges, out_fraction=%.3g",
                Wn.shape[0], Wn.nnz, out_fraction)
    return _transpose_and_rowsum(Wn)


def build_air_coupling(
    cfg: Config, top_n: int | None, n_nodes: int, out_fraction: float
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Load the cross-border air network for the chosen top-N airports as a coupling layer.

    Args:
        cfg: Project configuration.
        top_n: Airport cut-off; ``None`` uses the canonical ``global_air_network.parquet``,
            otherwise ``global_air_network_top<N>.parquet`` (produced by step 06 with the
            matching ``--top-n``).
        n_nodes: Total number of global nodes (matrix dimension).
        out_fraction: International FOI export fraction for gateway nodes.

    Returns:
        Tuple ``(WT, rowsum)`` for the air layer.

    Raises:
        FileNotFoundError: If the air-network file for the requested cut-off is missing.
    """
    name = "global_air_network.parquet" if top_n is None else f"global_air_network_top{top_n}.parquet"
    path = cfg.networks_dir / name
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found -- run scripts/06_global_air_network.py --top-n {top_n} first"
        )
    e = pd.read_parquet(path)
    W = sp.coo_matrix(
        (e["weight"].to_numpy(),
         (e["src_global_nodeid"].to_numpy(), e["dst_global_nodeid"].to_numpy())),
        shape=(n_nodes, n_nodes),
    ).tocsr()
    Wn = row_normalize(W, out_fraction)
    n_gateways = int((np.asarray(Wn.sum(axis=1)).ravel() > 0).sum())
    logger.info("abm.networks: air coupling (top-%s) %d edges across %d gateway nodes, out_fraction=%.3g",
                top_n, Wn.nnz, n_gateways, out_fraction)
    return _transpose_and_rowsum(Wn)
