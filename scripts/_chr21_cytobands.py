"""Standard UCSC hg38 cytoband stain colors.

Originally an external dependency; now vendored into the project.
"""

from __future__ import annotations

STAIN_COLORS: dict[str, str] = {
    "acen": "#D32F2F",
    "gpos": "#212121",
    "gpos100": "#212121",
    "gpos75": "#616161",
    "gpos66": "#757575",
    "gpos50": "#9E9E9E",
    "gpos33": "#BDBDBD",
    "gpos25": "#E0E0E0",
    "gneg": "#F5F5F5",
    "gvar": "#455A64",
    "stalk": "#90CAF9",
}


def acen_band_fill() -> str:
    """Fill color for centromeric (acen) bands."""
    return STAIN_COLORS["acen"]
