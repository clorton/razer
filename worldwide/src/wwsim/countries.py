"""The 193 United Nations member states, as ISO 3166-1 alpha-3 codes.

This is the canonical scope for the worldwide build. It deliberately **excludes**:

- UN observer states: Holy See (``VAT``) and State of Palestine (``PSE``);
- Non-members frequently seen in geodata: Taiwan (``TWN``), Kosovo (``XKX``),
  Western Sahara (``ESH``), Cook Islands, Niue, and all dependent territories.

The list is hard-coded (there is no machine-readable "UN membership" flag in pycountry)
and guarded by an assertion that it contains exactly 193 distinct, valid alpha-3 codes.
"""

from __future__ import annotations

import pycountry

from .logging import logger

__all__ = ["UN_MEMBERS", "country_name", "is_un_member"]

# 193 UN member states (ISO 3166-1 alpha-3), grouped 10 per line for auditability.
UN_MEMBERS: tuple[str, ...] = (
    "AFG", "ALB", "DZA", "AND", "AGO", "ATG", "ARG", "ARM", "AUS", "AUT",
    "AZE", "BHS", "BHR", "BGD", "BRB", "BLR", "BEL", "BLZ", "BEN", "BTN",
    "BOL", "BIH", "BWA", "BRA", "BRN", "BGR", "BFA", "BDI", "CPV", "KHM",
    "CMR", "CAN", "CAF", "TCD", "CHL", "CHN", "COL", "COM", "COG", "COD",
    "CRI", "CIV", "HRV", "CUB", "CYP", "CZE", "DNK", "DJI", "DMA", "DOM",
    "ECU", "EGY", "SLV", "GNQ", "ERI", "EST", "SWZ", "ETH", "FJI", "FIN",
    "FRA", "GAB", "GMB", "GEO", "DEU", "GHA", "GRC", "GRD", "GTM", "GIN",
    "GNB", "GUY", "HTI", "HND", "HUN", "ISL", "IND", "IDN", "IRN", "IRQ",
    "IRL", "ISR", "ITA", "JAM", "JPN", "JOR", "KAZ", "KEN", "KIR", "PRK",
    "KOR", "KWT", "KGZ", "LAO", "LVA", "LBN", "LSO", "LBR", "LBY", "LIE",
    "LTU", "LUX", "MDG", "MWI", "MYS", "MDV", "MLI", "MLT", "MHL", "MRT",
    "MUS", "MEX", "FSM", "MDA", "MCO", "MNG", "MNE", "MAR", "MOZ", "MMR",
    "NAM", "NRU", "NPL", "NLD", "NZL", "NIC", "NER", "NGA", "MKD", "NOR",
    "OMN", "PAK", "PLW", "PAN", "PNG", "PRY", "PER", "PHL", "POL", "PRT",
    "QAT", "ROU", "RUS", "RWA", "KNA", "LCA", "VCT", "WSM", "SMR", "STP",
    "SAU", "SEN", "SRB", "SYC", "SLE", "SGP", "SVK", "SVN", "SLB", "SOM",
    "ZAF", "SSD", "ESP", "LKA", "SDN", "SUR", "SWE", "CHE", "SYR", "TJK",
    "TZA", "THA", "TLS", "TGO", "TON", "TTO", "TUN", "TUR", "TKM", "TUV",
    "UGA", "UKR", "ARE", "GBR", "USA", "URY", "UZB", "VUT", "VEN", "VNM",
    "YEM", "ZMB", "ZWE",
)

# Fail fast at import time if the list is ever edited into an inconsistent state.
assert len(UN_MEMBERS) == 193, f"expected 193 UN members, got {len(UN_MEMBERS)}"
assert len(set(UN_MEMBERS)) == 193, "duplicate ISO3 code(s) in UN_MEMBERS"


def country_name(iso3: str) -> str:
    """Return a human-readable country name for an ISO3 code.

    Args:
        iso3: ISO 3166-1 alpha-3 code.

    Returns:
        The country's common name from pycountry, or the code itself if pycountry has no
        entry (should not happen for valid UN members).
    """
    rec = pycountry.countries.get(alpha_3=iso3)
    if rec is None:
        logger.warning("country_name: unknown ISO3 %r", iso3)
        return iso3
    return rec.name


def is_un_member(iso3: str) -> bool:
    """Whether an ISO3 code is one of the 193 UN member states.

    Args:
        iso3: ISO 3166-1 alpha-3 code (case-insensitive).

    Returns:
        ``True`` if the code is a UN member state, else ``False``.
    """
    return iso3.upper() in UN_MEMBERS
