"""Integration tests driving the CLI end-to-end (primary seam)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from bilan_immuno_twin.bilan import PANEL


def _run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "bilan_immuno_twin.cli", *args],
        capture_output=True,
        text=True,
    )


@pytest.fixture()
def bilan_file(tmp_path: Path) -> Path:
    payload = {
        "populations": {
            "CD8": 0.15,
            "Th1": 0.04,
            "Th2": 0.06,
            "Th17": 0.03,
            "B": 0.12,
            "NK": 0.09,
            "Treg": 0.02,
            "Monocytes": 0.11,
        }
    }
    path = tmp_path / "bilan.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_cli_bilan_emits_populations_as_json(bilan_file: Path) -> None:
    result = _run_cli("bilan", str(bilan_file))
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert set(parsed["populations"]) == set(PANEL)
    assert parsed["populations"]["CD8"] == 0.15


def test_cli_bilan_missing_population_prints_error_and_fails(bilan_file: Path) -> None:
    payload = json.loads(bilan_file.read_text(encoding="utf-8"))
    del payload["populations"]["NK"]
    bilan_file.write_text(json.dumps(payload), encoding="utf-8")
    result = _run_cli("bilan", str(bilan_file))
    assert result.returncode == 2
    assert "NK" in result.stderr


def test_cli_bilan_missing_file_fails(bilan_file: Path) -> None:
    result = _run_cli("bilan", str(bilan_file.parent / "absent.json"))
    assert result.returncode == 2
    assert "not found" in result.stderr


def test_cli_provides_help() -> None:
    result = _run_cli("--help")
    assert result.returncode == 0
    assert "bilan" in result.stdout


def test_cli_simulate_emits_trajectories_as_json(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "28")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert set(parsed["populations"]) == set(PANEL)
    cd8 = parsed["populations"]["CD8"]
    assert len(cd8["times"]) == len(cd8["values"])
    assert len(cd8["times"]) >= 2
    assert cd8["times"][0] == 0.0
    assert abs(cd8["times"][-1] - 28.0) < 1e-6
    assert parsed["t_units"] == "days"


def test_cli_simulate_adaptive_flag_stops_at_steady_state(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "400", "--adaptive")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert parsed["stopped_early"] is True
    cd8 = parsed["populations"]["CD8"]
    assert cd8["times"][-1] < 400.0


def test_cli_simulate_adaptive_flag_false_does_not_stop(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "100")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert parsed["stopped_early"] is False
    cd8 = parsed["populations"]["CD8"]
    assert abs(cd8["times"][-1] - 100.0) < 1e-6


def test_cli_simulate_all_populations_finite(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "28")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    for pop in PANEL:
        assert all(isinstance(v, (int, float)) for v in parsed["populations"][pop]["values"])