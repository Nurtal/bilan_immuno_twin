"""Unit tests for bilan parsing and validation (T1 pure seam)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bilan_immuno_twin.bilan import BilanError, PANEL, parse_bilan


def _valid_bilan(**overrides) -> dict:
    base = {
        "CD8": 0.15,
        "Th1": 0.04,
        "Th2": 0.06,
        "Th17": 0.03,
        "B": 0.12,
        "NK": 0.09,
        "Treg": 0.02,
        "Monocytes": 0.11,
    }
    base.update(overrides)
    return base


def _write(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_parse_valid_bilan_returns_all_panel_values(tmp_path: Path) -> None:
    path = _write(tmp_path / "bilan.json", {"populations": _valid_bilan()})
    result = parse_bilan(str(path))
    assert set(result) == set(PANEL)
    assert result["CD8"] == 0.15
    assert result["Treg"] == 0.02


def test_panel_contains_separate_th_subsets_and_core_populations() -> None:
    assert {"CD8", "Th1", "Th2", "Th17", "B", "NK", "Treg", "Monocytes"} == set(PANEL)


def test_missing_population_raises_clear_error(tmp_path: Path) -> None:
    path = _write(tmp_path / "bilan.json", {"populations": _valid_bilan(Treg=None)})
    payload = {"populations": {k: v for k, v in _valid_bilan().items() if k != "Treg"}}
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(BilanError, match="Treg"):
        parse_bilan(str(path))


def test_unknown_population_raises_error(tmp_path: Path) -> None:
    path = _write(tmp_path / "bilan.json", {"populations": _valid_bilan(Something=0.1)})
    with pytest.raises(BilanError, match="Something"):
        parse_bilan(str(path))


def test_non_numeric_population_raises_error(tmp_path: Path) -> None:
    path = _write(tmp_path / "bilan.json", {"populations": _valid_bilan(B="many")})
    with pytest.raises(BilanError, match="B"):
        parse_bilan(str(path))


def test_negative_population_raises_error(tmp_path: Path) -> None:
    path = _write(tmp_path / "bilan.json", {"populations": _valid_bilan(NK=-0.1)})
    with pytest.raises(BilanError, match="non-negative"):
        parse_bilan(str(path))


def test_missing_file_raises_error(tmp_path: Path) -> None:
    with pytest.raises(BilanError, match="not found"):
        parse_bilan(str(tmp_path / "nope.json"))


def test_malformed_json_raises_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text("{not json", encoding="utf-8")
    with pytest.raises(BilanError, match="invalid JSON"):
        parse_bilan(str(path))


def test_missing_populations_key_raises_error(tmp_path: Path) -> None:
    path = _write(tmp_path / "bilan.json", {"foo": 1})
    with pytest.raises(BilanError, match="'populations'"):
        parse_bilan(str(path))
