"""Parsing and validation of a bilan immunologique (cytometry snapshot)."""

from __future__ import annotations

import json
from pathlib import Path

PANEL = [
    "CD8",       # CD8+ T
    "Th1",       # CD4+ Th1
    "Th2",       # CD4+ Th2
    "Th17",      # CD4+ Th17
    "B",         # B cells
    "NK",        # NK cells
    "Treg",      # regulatory T
    "Monocytes", # monocytes
]

TOTAL = "total_cells"


class BilanError(ValueError):
    """Raised when a bilan file is missing, malformed or does not match the panel."""


def _validate(populations: dict, source: str) -> dict:
    missing = [name for name in PANEL if name not in populations]
    if missing:
        raise BilanError(
            f"{source}: missing population(s) for the panel cytometrique: {', '.join(missing)}"
        )

    extra = [name for name in populations if name not in PANEL]
    if extra:
        raise BilanError(
            f"{source}: unknown population(s) not in the panel cytometrique: {', '.join(extra)}"
        )

    for name in PANEL:
        value = populations[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise BilanError(f"{source}: population '{name}' must be numeric, got {value!r}")
        if value < 0:
            raise BilanError(f"{source}: population '{name}' must be non-negative, got {value}")

    return {name: float(populations[name]) for name in PANEL}


def parse_bilan(path: str | Path) -> dict:
    """Load a bilan immunologique file and validate it against the panel cytometrique.

    The file is JSON with a required "populations" object mapping each population
    of the panel to a non-negative numeric value (fraction of total immune cells).
    Returns a dict of structured, validated population values.
    """
    path = Path(path)
    if not path.is_file():
        raise BilanError(f"bilan file not found: {path}")

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise BilanError(f"{path}: invalid JSON: {exc}") from exc

    if not isinstance(raw, dict) or "populations" not in raw:
        raise BilanError(f"{path}: expected a JSON object with a 'populations' key")

    populations = raw["populations"]
    if not isinstance(populations, dict):
        raise BilanError(f"{path}: 'populations' must be an object mapping population to value")

    return _validate(populations, source=str(path))
