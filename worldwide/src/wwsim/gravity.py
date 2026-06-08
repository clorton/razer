"""Per-country intra-country gravity migration matrices over admin-2 nodes.

For each country we build a dense N×N matrix of connection strengths between its admin-2
units using the classic gravity model

    ``w_ij = k * P_i^a * P_j^b / D_ij^c``   (i != j; diagonal = 0)

where ``P`` is admin-2 population and ``D_ij`` is the great-circle distance (km) between
admin-2 centroids. These are the "per-country contagion migration matrices" -- one per
country -- that drive within-country spatial spread. The global cross-border coupling is
added separately by :mod:`wwsim.air_network`.

Matrices are written as ``.npz`` (the dense matrix + the ``global_nodeid`` index, so each
country's block lines up with the global node table) and as a Parquet edge list.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import Config, GravityParams
from .logging import logger

_EARTH_RADIUS_KM = 6371.0088


def haversine_matrix(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Pairwise great-circle distances (km) between points given in degrees.

    Args:
        lon: 1-D array of longitudes (degrees).
        lat: 1-D array of latitudes (degrees).

    Returns:
        Symmetric ``(N, N)`` matrix of great-circle distances in kilometres, zero on the
        diagonal.
    """
    lon_r = np.radians(np.asarray(lon, dtype=float))
    lat_r = np.radians(np.asarray(lat, dtype=float))
    dlon = lon_r[:, None] - lon_r[None, :]
    dlat = lat_r[:, None] - lat_r[None, :]
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat_r)[:, None] * np.cos(lat_r)[None, :] * np.sin(dlon / 2.0) ** 2
    a = np.clip(a, 0.0, 1.0)
    return 2.0 * _EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def gravity_matrix(
    pop: np.ndarray, lon: np.ndarray, lat: np.ndarray, params: GravityParams
) -> np.ndarray:
    """Compute a dense gravity connection matrix for one set of nodes.

    Args:
        pop: 1-D array of node populations.
        lon: 1-D array of node-centroid longitudes (degrees).
        lat: 1-D array of node-centroid latitudes (degrees).
        params: Gravity parameters (``k, a, b, c, min_distance_km``).

    Returns:
        ``(N, N)`` float64 matrix of connection weights with a zero diagonal. Non-finite
        entries (e.g. from zero distance) are set to zero.
    """
    pop = np.asarray(pop, dtype=float)
    n = len(pop)
    if n == 0:
        return np.zeros((0, 0))

    dist = haversine_matrix(lon, lat)
    dist = np.maximum(dist, params.min_distance_km)  # distance floor avoids blow-up

    pi = pop[:, None] ** params.a
    pj = pop[None, :] ** params.b
    with np.errstate(divide="ignore", invalid="ignore"):
        w = params.k * pi * pj / (dist**params.c)
    np.fill_diagonal(w, 0.0)
    w[~np.isfinite(w)] = 0.0
    return w


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    """Row-normalize a matrix so each row sums to 1 (rows that sum to 0 stay 0).

    Useful when a model wants per-source out-migration *probabilities* rather than raw
    gravity weights.

    Args:
        matrix: A square connection matrix.

    Returns:
        The row-normalized matrix (same shape).
    """
    rowsum = matrix.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        normed = np.where(rowsum > 0, matrix / rowsum, 0.0)
    return normed


def build_country_matrix(
    nodes: pd.DataFrame, params: GravityParams
) -> tuple[np.ndarray, np.ndarray]:
    """Build one country's gravity matrix from its slice of the global node table.

    Args:
        nodes: Rows of the global node table for a single country (must have
            ``global_nodeid, population, lon, lat``).
        params: Gravity parameters.

    Returns:
        Tuple of (``global_nodeid`` array of length N, ``(N, N)`` gravity matrix). The id
        array maps each matrix row/column to its global node id.
    """
    ids = nodes["global_nodeid"].to_numpy()
    matrix = gravity_matrix(
        nodes["population"].to_numpy(), nodes["lon"].to_numpy(), nodes["lat"].to_numpy(), params
    )
    return ids, matrix


def edge_list_from_matrix(ids: np.ndarray, matrix: np.ndarray) -> pd.DataFrame:
    """Convert a dense matrix to a non-zero directed edge list keyed by global node id.

    Args:
        ids: Global node ids for the matrix rows/columns.
        matrix: Dense connection matrix.

    Returns:
        DataFrame ``[src_global_nodeid, dst_global_nodeid, weight]`` for all non-zero entries.
    """
    src, dst = np.nonzero(matrix)
    return pd.DataFrame(
        {
            "src_global_nodeid": ids[src],
            "dst_global_nodeid": ids[dst],
            "weight": matrix[src, dst],
        }
    )


def build_all_country_matrices(
    global_nodes: pd.DataFrame, cfg: Config, save: bool = True
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Build (and optionally save) a gravity matrix for every country in the node table.

    Args:
        global_nodes: The global node table (from :mod:`wwsim.nodes`).
        cfg: Project configuration (gravity params + output dir).
        save: If ``True``, write ``<ISO3>_gravity.npz`` and ``<ISO3>_gravity.parquet`` per
            country under ``output/networks/``.

    Returns:
        Mapping ``iso3 -> (global_nodeid array, dense matrix)``.
    """
    cfg.networks_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for iso3, grp in global_nodes.groupby("iso3"):
        ids, matrix = build_country_matrix(grp, cfg.gravity)
        results[iso3] = (ids, matrix)
        if save:
            npz_path = cfg.networks_dir / f"{iso3}_gravity.npz"
            np.savez_compressed(npz_path, ids=ids, matrix=matrix.astype(np.float32))
            edge_list_from_matrix(ids, matrix).to_parquet(
                cfg.networks_dir / f"{iso3}_gravity.parquet", index=False
            )
    logger.info("gravity: built matrices for %d countries", len(results))
    return results


def country_matrices_to_global_sparse(
    matrices: dict[str, tuple[np.ndarray, np.ndarray]], n_nodes: int
) -> sp.csr_matrix:
    """Assemble all per-country matrices into one block-diagonal global sparse matrix.

    Args:
        matrices: Mapping ``iso3 -> (ids, dense matrix)`` from
            :func:`build_all_country_matrices`.
        n_nodes: Total number of global nodes (matrix dimension).

    Returns:
        An ``(n_nodes, n_nodes)`` CSR matrix containing each country's gravity block placed
        at its global-node indices (off-country entries are zero -- intra-country only).
    """
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    for ids, matrix in matrices.values():
        r, c = np.nonzero(matrix)
        rows.append(ids[r])
        cols.append(ids[c])
        vals.append(matrix[r, c])
    if not rows:
        return sp.csr_matrix((n_nodes, n_nodes))
    return sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_nodes, n_nodes),
    ).tocsr()
