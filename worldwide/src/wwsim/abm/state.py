"""Agent disease states, as small integer constants for a ``uint8`` state column.

Plain module-level ints (not an Enum) so numba njit kernels can close over them as
compile-time constants. SEIR uses S/E/I/R; an optional waning path (SEIRS) sends R back to
S. There is no DECEASED state -- this is a closed-population model (no births/deaths).
"""

from __future__ import annotations

SUSCEPTIBLE = 0
EXPOSED = 1
INFECTIOUS = 2
RECOVERED = 3

STATE_NAMES = {SUSCEPTIBLE: "S", EXPOSED: "E", INFECTIOUS: "I", RECOVERED: "R"}
