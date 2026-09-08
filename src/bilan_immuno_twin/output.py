"""CSV serialization of CLI outputs for downstream tooling (US #21, #22).

The JSON outputs carry the full structure; the CSV forms are flat tables that
drop straight into spreadsheets or pandas:

- a bilan -> one row per population (`population,value`)
- a simulated trajectory -> long format (`time,scenario,population,value`,
  plus `ci_lo,ci_hi` when confidence bands are available)
- a comparison -> `population,fold_change,direction`
- a score de réponse -> a single `therapy,score,...` row
- a calibration -> one row per population (`population,growth,ci_lo,ci_hi`)

When a command produces several tables (e.g. `simulate --perturb` emits
trajectories + comparison + score), they are concatenated separated by a blank
line, each beginning with its header.
"""

from __future__ import annotations

import csv
import io
from typing import Mapping, Sequence

from bilan_immuno_twin.bilan import PANEL


def _render(rows: Sequence[Sequence[str]]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerows(rows)
    return buffer.getvalue()


def populations_csv(populations: Mapping[str, float]) -> str:
    """Flat table of measured populations: ``population,value``."""
    return _render([["population", "value"], *[[pop, str(populations[pop])] for pop in PANEL]])


def trajectory_csv(
    trajectory: dict, scenario: str = "baseline", include_header: bool = True
) -> str:
    """Long-format trajectory table from a simulate output.

    ``trajectory`` is the ``{"populations": {pop: {"times", "values", ...}}}``
    shape emitted by the CLI. The ``ci_lo``/``ci_hi`` columns are added only
    when those fields are present (i.e. ``--with-ci``). Set ``include_header``
    to False to append another scenario to an already-headed table.
    """
    populations = trajectory["populations"]
    has_ci = "ci_lo" in populations[PANEL[0]] and "ci_hi" in populations[PANEL[0]]

    header = ["time", "scenario", "population", "value"]
    if has_ci:
        header += ["ci_lo", "ci_hi"]

    rows = [header] if include_header else []
    for pop in PANEL:
        entry = populations[pop]
        for i, value in enumerate(entry["values"]):
            row = [str(entry["times"][i]), scenario, pop, str(value)]
            if has_ci:
                row += [str(entry["ci_lo"][i]), str(entry["ci_hi"][i])]
            rows.append(row)
    return _render(rows)


def comparison_csv(comparison: Mapping[str, dict]) -> str:
    """What-if comparison table: ``population,fold_change,direction``."""
    rows = [["population", "fold_change", "direction"]]
    for pop in PANEL:
        fold = comparison[pop]["fold_change"]
        rows.append([pop, "" if fold is None else str(fold), comparison[pop]["direction"]])
    return _render(rows)


def response_csv(response: Mapping[str, object]) -> str:
    """Score de réponse row: a single ``therapy,score,...`` row.

    ``response`` is the response block emitted by ``simulate --perturb``.
    """
    header = [
        "therapy", "score", "interpretation", "reference_score",
        "is_differential", "unexpected_populations",
    ]
    row = [
        str(response["therapy"]),
        str(response["score"]),
        str(response["interpretation"]),
        str(response["reference_score"]),
        str(response["is_differential"]),
        ";".join(response["unexpected_populations"]),  # type: ignore[arg-type]
    ]
    return _render([header, row])


def calibration_csv(growth: Mapping[str, float], intervals: Mapping[str, Sequence[float]]) -> str:
    """Calibration table: ``population,growth,ci_lo,ci_hi``."""
    rows = [["population", "growth", "ci_lo", "ci_hi"]]
    for pop in PANEL:
        lo, hi = intervals[pop]
        rows.append([pop, str(growth[pop]), str(lo), str(hi)])
    return _render(rows)


def metadata_csv(meta: Mapping[str, object]) -> str:
    """Two-column table of key/value metadata (e.g. perturbation identity,
    calibration quality), where values may be numbers, strings or lists of
    strings flattened with ``;``."""
    rows = [["key", "value"]]
    for key, value in meta.items():
        if isinstance(value, str):
            text = value
        elif isinstance(value, (list, tuple)):
            text = ";".join(str(item) for item in value)
        else:
            text = str(value)
        rows.append([key, text])
    return _render(rows)