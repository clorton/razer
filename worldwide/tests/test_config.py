"""Tests for configuration loading.

Failure here means scripts run with wrong parameters or write to wrong paths.
"""

from __future__ import annotations

import pytest

from wwsim.config import Config, load_config


def test_defaults():
    """Given no overrides, when a Config is built, then 2015/admin-2 defaults hold."""
    c = Config()
    assert c.year == 2015
    assert c.adm_level == 2
    assert c.gravity.c == 2.0
    assert c.air.intra_country_air is False


def test_yaml_overlay(tmp_path):
    """Given a YAML overlay, when loaded, then only the named keys are overridden.

    A wrong overlay merge would silently ignore user parameters or clobber defaults.
    """
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text(
        "year: 2020\n"
        "gravity:\n  c: 1.5\n"
        "air:\n  top_n_airports: 250\n"
    )
    c = load_config(cfg_file)
    assert c.year == 2020
    assert c.gravity.c == 1.5
    assert c.gravity.k == 500.0  # untouched default preserved
    assert c.air.top_n_airports == 250


def test_missing_config_raises(tmp_path):
    """Given a non-existent path, when loaded, then FileNotFoundError is raised."""
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "nope.yaml")
