"""Worldwide agent-based SEIR model built on laser-core, with a custom sparse FOI.

- :class:`~wwsim.abm.model.WorldSEIR` / :class:`~wwsim.abm.model.SEIRParams` -- the model.
- :mod:`wwsim.abm.networks` -- sparse intra-country (+ optional air) coupling layers.
- :mod:`wwsim.abm.components` -- single-``uint8``-timer progression and the custom FOI.
"""

from __future__ import annotations

from .model import SEIRParams, WorldSEIR

__all__ = ["SEIRParams", "WorldSEIR"]
