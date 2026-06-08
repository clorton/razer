"""Per-tick ABM components: single-timer progression and custom sparse-network FOI.

Two components run each tick, in this order (the razer ordering -- progression *before*
transmission -- so a just-infected agent is not decremented the same tick and every state
realizes its full duration, giving ``R0 = beta * infectious_period`` with no off-by-one):

1. :class:`Progression` -- one numba pass over agents that **branches on the entry state**
   and reuses a single ``uint8`` ``timer`` for E→I→(R), plus an optional R→S waning. Because
   it branches on the state each agent had at the top of the tick, an agent that turns E→I
   here is not also stepped as I, so one timer serves all timed states.
2. :class:`Transmission` -- the **custom force of infection**. It mirrors laser-generic's
   ``ft = beta*I/N`` then network coupling, but does the coupling with **sparse** matrix-
   vector products instead of a dense ``N×N`` multiply: a node exports a fraction of its FOI
   to its neighbours (intra-country gravity) and, optionally, to gateway nodes abroad (air).

Node compartment counts (`nodes.S/E/I/R[tick]`) are carried forward `t -> t+1` by the model
and updated here by the per-node transition deltas the kernels accumulate.
"""

from __future__ import annotations

import numba as nb
import numpy as np
import scipy.sparse as sp

from .state import EXPOSED, INFECTIOUS, RECOVERED, SUSCEPTIBLE


@nb.njit(nogil=True, parallel=True, cache=True)
def nb_progress(state, timer, nodeid, inf_dur, wan_dur, has_waning,
                to_infectious, to_recovered, to_susceptible):
    """Advance timed states by one day, branching on each agent's entry state.

    Args:
        state: uint8 per-agent state array (mutated).
        timer: uint8 per-agent countdown (mutated; reused across E/I/R).
        nodeid: uint16 per-agent node id.
        inf_dur: infectious-period length (days) set when an agent turns infectious.
        wan_dur: waning-immunity length (days) set on recovery when ``has_waning``.
        has_waning: whether R wanes back to S (SEIRS) or is terminal (SEIR).
        to_infectious / to_recovered / to_susceptible: ``(threads, nodes)`` scratch arrays
            into which per-node transition counts are accumulated (summed by the caller).
    """
    for i in nb.prange(len(state)):
        s = state[i]
        tid = nb.get_thread_id()
        nid = nodeid[i]
        if s == EXPOSED:
            timer[i] -= 1
            if timer[i] == 0:
                state[i] = INFECTIOUS
                timer[i] = inf_dur
                to_infectious[tid, nid] += 1
        elif s == INFECTIOUS:
            timer[i] -= 1
            if timer[i] == 0:
                state[i] = RECOVERED
                timer[i] = wan_dur if has_waning else 0
                to_recovered[tid, nid] += 1
        elif has_waning and s == RECOVERED:
            timer[i] -= 1
            if timer[i] == 0:
                state[i] = SUSCEPTIBLE
                to_susceptible[tid, nid] += 1
    return


@nb.njit(nogil=True, parallel=True, cache=True)
def nb_infect(state, nodeid, prob, inc_dur, timer, newly_infected):
    """Expose susceptible agents with per-node probability ``prob[node]``.

    Args:
        state: uint8 per-agent state (mutated S→E).
        nodeid: uint16 per-agent node id.
        prob: per-node infection probability for this tick.
        inc_dur: incubation length (days) set on the newly exposed.
        timer: uint8 per-agent countdown (mutated).
        newly_infected: ``(threads, nodes)`` scratch array of new exposures per node.
    """
    for i in nb.prange(len(state)):
        if state[i] == SUSCEPTIBLE:
            nid = nodeid[i]
            if np.random.random() < prob[nid]:
                state[i] = EXPOSED
                timer[i] = inc_dur
                newly_infected[nb.get_thread_id(), nid] += 1
    return


class Progression:
    """E→I→R (+ optional R→S) using one reused ``uint8`` timer."""

    def __init__(self, model):
        """Bind to the model and read durations from its params.

        Args:
            model: The :class:`wwsim.abm.model.WorldSEIR` instance.
        """
        self.model = model
        p = model.params
        self.inf_dur = int(p.infectious_period)
        self.wan_dur = int(getattr(p, "waning_period", 0) or 0)
        self.has_waning = self.wan_dur > 0

    def step(self, tick: int) -> None:
        """Run the progression pass and apply per-node E/I/R deltas to ``tick+1``."""
        m = self.model
        threads = nb.get_num_threads()
        nn = m.nodes.count
        to_inf = np.zeros((threads, nn), dtype=np.int32)
        to_rec = np.zeros((threads, nn), dtype=np.int32)
        to_sus = np.zeros((threads, nn), dtype=np.int32)

        nb_progress(m.people.state, m.people.timer, m.people.nodeid,
                    self.inf_dur, self.wan_dur, self.has_waning, to_inf, to_rec, to_sus)

        ni = to_inf.sum(axis=0)   # E -> I
        nr = to_rec.sum(axis=0)   # I -> R
        ns = to_sus.sum(axis=0)   # R -> S (waning)

        m.nodes.E[tick + 1] -= ni
        m.nodes.I[tick + 1] += ni
        m.nodes.I[tick + 1] -= nr
        m.nodes.R[tick + 1] += nr
        if self.has_waning:
            m.nodes.R[tick + 1] -= ns
            m.nodes.S[tick + 1] += ns
        m.nodes.newly_infectious[tick] = ni
        m.nodes.newly_recovered[tick] = nr
        return


class Transmission:
    """Custom force of infection: ``beta*I/N`` plus sparse intra (+ optional air) coupling.

    The coupling reproduces laser-generic's dense ``ft + Wᵀft − ft·rowsum`` with sparse
    matrix-vector products. The network layers are pre-transposed (``WT``) with cached row
    sums by :mod:`wwsim.abm.networks`.
    """

    def __init__(self, model, intra, air=None):
        """Bind to the model and store the coupling layers.

        Args:
            model: The :class:`wwsim.abm.model.WorldSEIR` instance.
            intra: ``(WT, rowsum)`` for the intra-country block-diagonal coupling.
            air: Optional ``(WT, rowsum)`` for the cross-border air coupling.
        """
        self.model = model
        self.beta = float(model.params.beta)
        self.inc_dur = int(model.params.incubation_period)
        self.intra_WT, self.intra_rowsum = intra
        self.air_WT, self.air_rowsum = air if air is not None else (None, None)
        # Optional per-country, per-tick beta multiplier (e.g. empirical Rt forcing) and the
        # node->country index that expands it to per-node. None => constant beta.
        self.forcing = getattr(model, "forcing", None)
        self.node_country_idx = getattr(model, "node_country_idx", None)

    def step(self, tick: int) -> None:
        """Compute per-node FOI, apply network coupling, and expose susceptibles."""
        m = self.model
        # Infectious census after this tick's progression (E-entry: new infectious counted,
        # recoveries excluded) -> contributes on exactly `infectious_period` ticks => beta*D.
        infectious = m.nodes.I[tick + 1].astype(np.float64)
        # Per-node transmission rate: beta, optionally scaled by this country's Rt multiplier
        # for this tick (each node looks up its country's factor).
        if self.forcing is not None:
            beta_node = self.beta * self.forcing[tick][self.node_country_idx]
            ft0 = beta_node * infectious / m.node_pop
        else:
            ft0 = self.beta * infectious / m.node_pop  # node_pop is clamped >= 1

        # Network coupling: incoming (Wᵀ·ft) minus outgoing (ft·rowsum). Sparse, O(nnz).
        ft = ft0 + self.intra_WT.dot(ft0) - ft0 * self.intra_rowsum
        if self.air_WT is not None:
            ft = ft + self.air_WT.dot(ft0) - ft0 * self.air_rowsum

        m.nodes.forces[tick] = ft.astype(np.float32)
        prob = -np.expm1(-ft)  # FOI -> per-node probability of infection this tick

        threads = nb.get_num_threads()
        newly = np.zeros((threads, m.nodes.count), dtype=np.int32)
        nb_infect(m.people.state, m.people.nodeid, prob, self.inc_dur, m.people.timer, newly)
        ni = newly.sum(axis=0)

        m.nodes.S[tick + 1] -= ni
        m.nodes.E[tick + 1] += ni
        m.nodes.newly_infected[tick] = ni
        return


def coupling_is_sparse(WT) -> bool:
    """Whether a coupling layer is a sparse matrix (helper for sanity checks/tests)."""
    return sp.issparse(WT)
