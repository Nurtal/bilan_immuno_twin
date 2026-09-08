"""CSV serialization of CLI outputs for downstream tooling (US #21, #22)."""

from __future__ import annotations

from bilan_immuno_twin.bilan import PANEL
from bilan_immuno_twin.output import (
    calibration_csv,
    comparison_csv,
    metadata_csv,
    populations_csv,
    response_csv,
    trajectory_csv,
)


def test_populations_csv_has_header_then_one_row_per_population() -> None:
    populations = {pop: float(i) / 100 for i, pop in enumerate(PANEL)}
    body = populations_csv(populations)
    lines = body.strip().splitlines()
    assert lines[0] == "population,value"
    assert len(lines) == len(PANEL) + 1
    assert lines[1] == "CD8,0.0"
    assert "Monocytes," in lines[-1]


def test_trajectory_csv_is_long_format_covering_every_step_and_population() -> None:
    times = [0.0, 1.0, 2.0]
    values = {pop: [0.1 * i for i in range(3)] for pop in PANEL}
    body = trajectory_csv({"populations": {
        pop: {"times": times, "values": values[pop]} for pop in PANEL
    }}, scenario="baseline")
    lines = body.strip().splitlines()
    assert lines[0] == "time,scenario,population,value"
    assert len(lines) == len(PANEL) * len(times) + 1
    assert lines[1] == "0.0,baseline,CD8,0.0"
    assert lines[2] == "1.0,baseline,CD8,0.1"


def test_trajectory_csv_carries_the_scenario_label() -> None:
    body = trajectory_csv({"populations": {
        pop: {"times": [0.0], "values": [0.0]} for pop in PANEL
    }}, scenario="perturbed")
    assert "0.0,perturbed,CD8,0.0" in body


def test_trajectory_csv_includes_ci_columns_only_when_available() -> None:
    plain = {"populations": {pop: {"times": [0.0], "values": [0.0]} for pop in PANEL}}
    with_ci = {"populations": {pop: {
        "times": [0.0], "values": [0.0], "ci_lo": [0.0], "ci_hi": [0.0],
    } for pop in PANEL}}
    assert trajectory_csv(plain).splitlines()[0] == "time,scenario,population,value"
    rich = trajectory_csv(with_ci).splitlines()[0]
    assert rich == "time,scenario,population,value,ci_lo,ci_hi"


def test_calibration_csv_has_estimate_and_interval_per_population() -> None:
    growth = {pop: 0.3 for pop in PANEL}
    intervals = {pop: (0.2, 0.4) for pop in PANEL}
    body = calibration_csv(growth, intervals)
    lines = body.strip().splitlines()
    assert lines[0] == "population,growth,ci_lo,ci_hi"
    assert len(lines) == len(PANEL) + 1
    assert lines[1] == "CD8,0.3,0.2,0.4"


def test_calibration_csv_accepts_list_bands_directly() -> None:
    growth = {pop: 0.3 for pop in PANEL}
    bands = {pop: [0.2, 0.4] for pop in PANEL}  # same shape as the JSON output
    assert calibration_csv(growth, bands).splitlines()[1] == "CD8,0.3,0.2,0.4"


def test_comparison_csv_one_row_per_population_with_direction() -> None:
    comparison = {
        pop: {"fold_change": 1.2 if pop != "B" else None, "direction": "up" if pop != "B" else "flat"}
        for pop in PANEL
    }
    body = comparison_csv(comparison)
    lines = body.strip().splitlines()
    assert lines[0] == "population,fold_change,direction"
    assert len(lines) == len(PANEL) + 1
    assert lines[1] == "CD8,1.2,up"
    assert "B,,flat" in lines  # no baseline means an empty cell, not a crash


def test_response_csv_emits_score_and_interpretation() -> None:
    response = {
        "therapy": "anti-PD1",
        "score": 1.5,
        "interpretation": "favorable",
        "reference_score": 0.9,
        "is_differential": True,
        "unexpected_populations": ["B", "Th2"],
        "folds": {},
    }
    body = response_csv(response)
    lines = body.strip().splitlines()
    assert lines[0] == (
        "therapy,score,interpretation,reference_score,is_differential,unexpected_populations"
    )
    assert lines[1] == "anti-PD1,1.5,favorable,0.9,True,B;Th2"


def test_metadata_csv_is_two_column_key_value() -> None:
    body = metadata_csv({"n_bootstrap": 100, "converged": True, "mse": 0.001})
    lines = body.strip().splitlines()
    assert lines[0] == "key,value"
    assert lines[1:4] == ["n_bootstrap,100", "converged,True", "mse,0.001"]