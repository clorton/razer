"""The worldwide agent-based SEIR model.

One agent per (subsampled) person, carrying just three columns -- ``state`` (uint8),
``nodeid`` (uint16), ``timer`` (uint8) = 4 bytes/agent -- so full resolution (~7.3 B agents
≈ 29 GB) fits in >32 GB RAM, and a ``subsample`` divisor makes test runs light.

The model holds a `people` LaserFrame (agents) and a `nodes` LaserFrame whose S/E/I/R are
``(nticks+1, n_nodes)`` count arrays carried forward each tick. It is deliberately lean: it
does **not** build laser-generic's dense gravity network (8 GB at this scale) -- the sparse
coupling layers are supplied by :mod:`wwsim.abm.networks` and consumed by the custom
:class:`wwsim.abm.components.Transmission`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from laser.core import LaserFrame
from laser.core.random import seed as set_seed

from ..logging import logger
from .components import Progression, Transmission
from .state import INFECTIOUS, SUSCEPTIBLE


@dataclass
class SEIRParams:
    """Parameters for a worldwide SEIR run.

    Attributes:
        nticks: Number of daily time steps.
        beta: Transmission rate per day (``R0 = beta * infectious_period`` in a naive,
            well-mixed node).
        incubation_period: Days from exposure to infectiousness (E→I), ``1..255``.
        infectious_period: Days infectious before recovery (I→R), ``1..255``.
        waning_period: Days of immunity before R→S; ``0`` = no waning (plain SEIR).
        subsample: Divide each node's population by this factor (``8`` ≈ 1/8 scale).
        intra_out_fraction: Fraction of a node's FOI exported to its country's other nodes.
        use_air: Whether to add the cross-border air coupling.
        air_top_n: Airport cut-off for the air layer (matches a built ``*_top<N>`` network);
            ``None`` uses the canonical air network.
        air_out_fraction: International FOI export fraction for gateway nodes.
        seed_count: Number of agents to seed infectious at start.
        seed_nodeid: Global node id to seed; ``None`` = the most populous node.
        prng_seed: RNG seed for reproducibility.
    """

    nticks: int = 180
    beta: float = 0.3
    incubation_period: int = 4
    infectious_period: int = 6
    waning_period: int = 0
    subsample: int = 1
    intra_out_fraction: float = 0.1
    use_air: bool = False
    air_top_n: int | None = None
    air_out_fraction: float = 0.02
    seed_count: int = 100
    seed_nodeid: int | None = None
    prng_seed: int = 20260101


class WorldSEIR:
    """Worldwide agent-based SEIR over admin-2 nodes with sparse network coupling."""

    def __init__(self, nodes_df: pd.DataFrame, params: SEIRParams, intra, air=None, forcing=None):
        """Build agents and node-count buffers, seed infections, and wire components.

        Args:
            nodes_df: Global node table (must be sorted so row ``k`` is ``global_nodeid==k``;
                needs ``population`` and ``iso3``). Defines node order = network index.
            params: Run parameters.
            intra: ``(WT, rowsum)`` intra-country coupling from :mod:`wwsim.abm.networks`.
            air: Optional ``(WT, rowsum)`` air coupling.
            forcing: Optional ``(m, iso3_order)`` where ``m`` is an ``(nticks, n_countries)``
                per-country, per-tick transmission multiplier (e.g. an Rt forcing from
                :func:`wwsim.covid_reference.rt_forcing`) and ``iso3_order`` is its column
                order. Each admin-2 node uses its country's factor.

        Raises:
            ValueError: If ``nodes_df`` is not indexed 0..N-1 by ``global_nodeid``, or the
                forcing's ``iso3_order`` does not cover every node's country.
        """
        self.params = params
        self.scenario = nodes_df.reset_index(drop=True)
        if "global_nodeid" in self.scenario.columns and not np.array_equal(
            self.scenario["global_nodeid"].to_numpy(), np.arange(len(self.scenario))
        ):
            raise ValueError("nodes_df must be sorted so row k == global_nodeid k")

        set_seed(params.prng_seed)
        np.random.seed(params.prng_seed)

        n_nodes = len(self.scenario)
        pops = np.maximum(np.round(self.scenario["population"].to_numpy() / params.subsample), 0).astype(np.int64)
        num_agents = int(pops.sum())
        logger.info("WorldSEIR: %d nodes, %d agents (subsample=%d, %.2f GB for 4B/agent)",
                    n_nodes, num_agents, params.subsample, 4 * num_agents / 1e9)

        # --- agents (one per subsampled person) ---
        self.people = LaserFrame(num_agents, num_agents)
        self.people.add_scalar_property("state", dtype=np.uint8, default=SUSCEPTIBLE)
        self.people.add_scalar_property("nodeid", dtype=np.uint16)
        self.people.add_scalar_property("timer", dtype=np.uint8, default=0)
        # Agents are laid out node-by-node, so node k owns the contiguous block
        # [offset[k] : offset[k+1]).
        offsets = np.zeros(n_nodes + 1, dtype=np.int64)
        np.cumsum(pops, out=offsets[1:])
        # Fill node ids in place. We deliberately avoid
        #   np.repeat(np.arange(n_nodes, dtype=np.int64), pops).astype(np.uint16)
        # because at low subsampling (billions of agents) that builds a ~29 GB int64
        # temporary (plus a 7 GB uint16 temporary) -- a ~3x peak-RAM spike at init that OOMs
        # in a container with a hard cgroup memory limit and little swap (macOS hides it via
        # memory compression + swap). Writing block-by-block into the preallocated uint16
        # column uses no large temporary.
        nid = self.people.nodeid
        for k in range(n_nodes):
            nid[offsets[k] : offsets[k + 1]] = k

        # --- node compartment buffers, carried forward each tick ---
        self.nodes = LaserFrame(n_nodes)
        for st in ("S", "E", "I", "R"):
            self.nodes.add_vector_property(st, params.nticks + 1, dtype=np.int32)
        for flow in ("forces",):
            self.nodes.add_vector_property(flow, params.nticks + 1, dtype=np.float32)
        for flow in ("newly_infected", "newly_infectious", "newly_recovered"):
            self.nodes.add_vector_property(flow, params.nticks + 1, dtype=np.int32)
        self.nodes.S[0] = pops.astype(np.int32)

        # Denominator N per node (constant: agents never move between nodes), clamped >= 1.
        self.node_pop = np.maximum(pops.astype(np.float64), 1.0)

        self._seed(params, pops, offsets)

        # Optional per-country, per-tick transmission multiplier, expanded to a per-node
        # lookup so each admin-2 node applies its country's factor in the FOI.
        self.forcing = None
        self.node_country_idx = None
        if forcing is not None:
            m_ct, iso3_order = forcing
            order = {iso: k for k, iso in enumerate(iso3_order)}
            nci = self.scenario["iso3"].map(order)
            if nci.isna().any():
                raise ValueError("forcing iso3_order must cover every node's country")
            self.node_country_idx = nci.to_numpy().astype(np.int64)
            self.forcing = np.asarray(m_ct, dtype=np.float32)
            logger.info("WorldSEIR: Rt forcing enabled (%d ticks x %d countries)", *self.forcing.shape)

        self.components = [Progression(self), Transmission(self, intra, air)]
        return

    def _seed(self, params: SEIRParams, pops: np.ndarray, offsets: np.ndarray) -> None:
        """Seed initial infectious agents in one node and update its counts."""
        node = params.seed_nodeid if params.seed_nodeid is not None else int(np.argmax(pops))
        count = min(params.seed_count, int(pops[node]))
        if count <= 0:
            logger.warning("WorldSEIR: seed node %d has no agents; nothing seeded", node)
            return
        start = offsets[node]
        idx = np.arange(start, start + count)  # first `count` agents in the node's block
        self.people.state[idx] = INFECTIOUS
        self.people.timer[idx] = params.infectious_period
        self.nodes.S[0, node] -= count
        self.nodes.I[0, node] += count
        logger.info("WorldSEIR: seeded %d infectious in node %d (%s)", count, node,
                    self.scenario.get("adm2_name", pd.Series(["?"] * len(self.scenario))).iloc[node])
        return

    def _initialize_flows(self, tick: int) -> None:
        """Carry each compartment's counts forward ``t -> t+1`` before components apply deltas."""
        for st in ("S", "E", "I", "R"):
            arr = getattr(self.nodes, st)
            arr[tick + 1] = arr[tick]
        return

    def run(self) -> None:
        """Run all ticks: carry forward, then progression, then transmission."""
        for tick in range(self.params.nticks):
            self._initialize_flows(tick)
            for c in self.components:
                c.step(tick)
        logger.info("WorldSEIR: completed %d ticks", self.params.nticks)
        return

    def save_history(self, path, compartments=("S", "E", "I", "R")) -> None:
        """Save the full per-node, per-tick compartment counts for later animation/analysis.

        Writes a compressed ``.npz`` holding each requested compartment as an
        ``(nticks+1, n_nodes)`` int32 array (columns aligned to ``global_nodeid`` 0..N-1),
        plus ``node_pop`` (the subsampled per-node denominator) and ``nticks``. Early ticks
        are mostly zeros, so compression keeps the file modest.

        Args:
            path: Output ``.npz`` path.
            compartments: Which compartments to store (default all of S/E/I/R).
        """
        arrays = {st: getattr(self.nodes, st) for st in compartments}
        np.savez_compressed(path, node_pop=self.node_pop, nticks=np.int64(self.params.nticks), **arrays)
        logger.info("WorldSEIR: saved per-node history (%s) -> %s",
                    ",".join(compartments), path)
        return

    def totals(self) -> pd.DataFrame:
        """Global time series of compartment totals (summed over all nodes).

        Returns:
            DataFrame indexed by tick with columns ``S, E, I, R`` and ``incidence`` (daily
            new exposures summed over nodes).
        """
        df = pd.DataFrame({st: getattr(self.nodes, st).sum(axis=1) for st in ("S", "E", "I", "R")})
        df["incidence"] = self.nodes.newly_infected.sum(axis=1)
        df.index.name = "tick"
        return df
