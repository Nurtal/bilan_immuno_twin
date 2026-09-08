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