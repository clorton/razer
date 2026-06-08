"""Combine per-mode networks into one global migration matrix (pluggable, rail-ready).

Every transport mode is reduced to the same currency: a directed edge list over global
admin-2 node ids, ``[src_global_nodeid, dst_global_nodeid, weight]``. A :class:`ModeNetwork`
wraps such an edge list with a name and a scale factor. :func:`combine_modes` sums them
into one sparse adjacency matrix.

This is the seam the goal asks for: today we combine

- **intra-country gravity** (block-diagonal; :mod:`wwsim.gravity`) and
- **inter-country air** (cross-border; :mod:`wwsim.air_network`),

and tomorrow a **rail** mode (:mod:`wwsim.rail`) or a licensed **OAG** air feed
(:mod:`wwsim.oag`) is just another :class:`ModeNetwork` added to the list -- no downstream
change. Per-mode ``scale`` makes otherwise-incommensurate units (gravity weights vs.
seats/yr) comparable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import Config
from .gravity import edge_list_from_matrix
from .logging import logger


@dataclass
class ModeNetwork:
    """One transport mode's contribution to the global network.

    Attributes:
        name: Mode label (e.g. ``"gravity"``, ``"air"``, ``"rail"``).
        edges: Directed edge list with ``src_global_nodeid, dst_global_nodeid, weight``.
        scale: Multiplier applied to this mode's weights when combining (default 1.0).
    """

    name: str
    edges: pd.DataFrame
    scale: float = 1.0


def combine_modes(modes: list[ModeNetwork], n_nodes: int) -> sp.csr_matrix:
    """Sum several modes' edge lists into one global sparse adjacency matrix.

    Args:
        modes: List of :class:`ModeNetwork`. Edge lists may overlap; weights add.
        n_nodes: Total number of global nodes (matrix dimension).

    Returns:
        An ``(n_nodes, n_nodes)`` CSR matrix; entry ``[i, j]`` is the scaled sum over modes
        of the i->j weight.
    """
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    for mode in modes:
        if len(mode.edges) == 0:
            continue
        rows.append(mode.edges["src_global_nodeid"].to_numpy())
        cols.append(mode.edges["dst_global_nodeid"].to_numpy())
        vals.append(mode.edges["weight"].to_numpy() * mode.scale)
        logger.info("combine_modes: %s contributes %d edges (scale=%g)",
                    mode.name, len(mode.edges), mode.scale)
    if not rows:
        return sp.csr_matrix((n_nodes, n_nodes))
    return sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_nodes, n_nodes),
    ).tocsr()


def gravity_edges_from_matrices(
    matrices: dict[str, tuple[np.ndarray, np.ndarray]]
) -> pd.DataFrame:
    """Flatten per-country gravity matrices into a single global edge list.

    Args:
        matrices: Mapping ``iso3 -> (ids, dense matrix)`` from
            :func:`wwsim.gravity.build_all_country_matrices`.

    Returns:
        Concatenated edge list ``[src_global_nodeid, dst_global_nodeid, weight]``.
    """
    frames = [edge_list_from_matrix(ids, m) for ids, m in matrices.values()]
    if not frames:
        return pd.DataFrame(columns=["src_global_nodeid", "dst_global_nodeid", "weight"])
    return pd.concat(frames, ignore_index=True)


def build_combined_network(
    n_nodes: int,
    gravity_matrices: dict[str, tuple[np.ndarray, np.ndarray]],
    air_edges: pd.DataFrame,
    cfg: Config,
    gravity_scale: float = 1.0,
    air_scale: float = 1.0,
    extra_modes: list[ModeNetwork] | None = None,
) -> sp.csr_matrix:
    """Assemble the combined global matrix from gravity + air (+ any extra modes).

    Args:
        n_nodes: Total number of global nodes.
        gravity_matrices: Per-country gravity matrices.
        air_edges: Global cross-border air edges (from :func:`wwsim.air_network.build_air_network`).
        cfg: Project configuration.
        gravity_scale: Scale for the intra-country gravity contribution.
        air_scale: Scale for the inter-country air contribution.
        extra_modes: Optional additional modes (e.g. rail) to add to the sum.

    Returns:
        The combined ``(n_nodes, n_nodes)`` CSR migration matrix.
    """
    modes = [
        ModeNetwork("gravity", gravity_edges_from_matrices(gravity_matrices), gravity_scale),
        ModeNetwork(
            "air",
            air_edges.rename(columns={}).loc[:, ["src_global_nodeid", "dst_global_nodeid", "weight"]],
            air_scale,
        ),
    ]
    if extra_modes:
        modes.extend(extra_modes)
    combined = combine_modes(modes, n_nodes)
    logger.info(
        "network: combined matrix %dx%d, %d nonzeros (intra block-diagonal + inter air)",
        n_nodes, n_nodes, combined.nnz,
    )
    return combined


def save_combined_network(
    matrix: sp.csr_matrix, cfg: Config, stem: str = "global_combined_network"
) -> None:
    """Persist the combined global matrix as a sparse ``.npz``.

    Args:
        matrix: Combined CSR matrix.
        cfg: Project configuration.
        stem: Output filename stem (default ``"global_combined_network"``); top-N variants
            use e.g. ``"global_combined_network_top250"`` so several can coexist.
    """
    cfg.networks_dir.mkdir(parents=True, exist_ok=True)
    sp.save_npz(cfg.networks_dir / f"{stem}.npz", matrix)
    logger.info("network: saved %s.npz", stem)
