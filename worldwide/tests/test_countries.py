"""Tests for the UN member-state list.

Failure here means the simulation scope is wrong -- a missing/extra country would silently
add or drop an entire nation's nodes and edges from the worldwide model.
"""

from __future__ import annotations

import pycountry

from wwsim.countries import UN_MEMBERS, country_name, is_un_member


def test_exactly_193_distinct_members():
    """Given the UN list, when counted, then there are exactly 193 distinct codes.

    The UN has 193 member states; duplicates or a wrong count corrupt the scope.
    """
    assert len(UN_MEMBERS) == 193
    assert len(set(UN_MEMBERS)) == 193


def test_all_codes_are_valid_iso3():
    """Given each code, when looked up in pycountry, then all resolve to a country."""
    invalid = [c for c in UN_MEMBERS if pycountry.countries.get(alpha_3=c) is None]
    assert invalid == []


def test_observers_and_nonmembers_excluded():
    """Given known non-members, when checked, then none are in the list.

    Holy See and Palestine are observers; Taiwan/Kosovo/Western Sahara are non-members.
    """
    for code in ("VAT", "PSE", "TWN", "XKX", "ESH"):
        assert code not in UN_MEMBERS


def test_is_un_member_is_case_insensitive():
    """Given a lowercase code, when is_un_member is called, then membership is detected."""
    assert is_un_member("nga") is True
    assert is_un_member("XKX") is False


def test_country_name_resolves_known_code():
    """Given a valid code, when country_name is called, then a human name is returned."""
    assert "Nigeria" in country_name("NGA")
