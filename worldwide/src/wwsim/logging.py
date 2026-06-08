"""Centralized logging for the ``wwsim`` package.

Every module imports ``logger`` from here (``from .logging import logger``) and logs
internal actions at ``INFO`` level so a long, mostly-unattended worldwide data build
leaves a readable trail.

The logger writes to ``stderr`` with a compact timestamped format. Importing this module
configures the handler exactly once (idempotent), so repeated imports across modules and
scripts never attach duplicate handlers.
"""

from __future__ import annotations

import logging
import sys

__all__ = ["logger", "get_logger"]

_LOGGER_NAME = "wwsim"
_FORMAT = "%(asctime)s %(levelname)-7s %(name)s: %(message)s"
_DATEFMT = "%H:%M:%S"


def _configure() -> logging.Logger:
    """Create and configure the package logger exactly once.

    Returns:
        The configured ``wwsim`` logger.
    """
    log = logging.getLogger(_LOGGER_NAME)
    # Only attach a handler the first time -- guards against duplicate lines when
    # several modules import this at different times.
    if not log.handlers:
        handler = logging.StreamHandler(stream=sys.stderr)
        handler.setFormatter(logging.Formatter(fmt=_FORMAT, datefmt=_DATEFMT))
        log.addHandler(handler)
        log.setLevel(logging.INFO)
        # Do not propagate to the root logger; we own this stream.
        log.propagate = False
    return log


logger = _configure()


def get_logger(suffix: str | None = None) -> logging.Logger:
    """Return the package logger, optionally a named child of it.

    Args:
        suffix: Optional child-logger name (e.g. ``"flights"``). Child loggers share the
            parent's handler and level, so output stays uniform.

    Returns:
        The ``wwsim`` logger, or ``wwsim.<suffix>`` if a suffix is given.
    """
    if suffix:
        return logger.getChild(suffix)
    return logger
