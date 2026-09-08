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


def test_cli_simulate_with_perturb_emits_comparison(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "28", "--perturb", "anti-PD1")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert parsed["perturbation"] == "anti-PD1"
    assert parsed["parameter_changes"]["growth:CD8"]["delta"] == 0.40
    assert set(parsed["unperturbed"]["populations"]) == set(PANEL)
    assert set(parsed["perturbed"]["populations"]) == set(PANEL)
    assert set(parsed["comparison"]) == set(PANEL)
    assert parsed["comparison"]["CD8"]["direction"] in ("up", "down", "flat")


def test_cli_simulate_unknown_perturb_fails_clearly(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--perturb", "vaccine-x")
    assert result.returncode == 2
    assert "unknown perturbation" in result.stderr


def test_cli_simulate_emits_response_score_and_flags(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "28", "--perturb", "anti-PD1")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    response = parsed["response"]
    assert response["therapy"] == "anti-PD1"
    assert isinstance(response["score"], float)
    assert isinstance(response["reference_score"], float)
    assert response["interpretation"] in ("favorable", "neutral", "unfavorable")
    assert isinstance(response["is_differential"], bool)
    assert isinstance(response["unexpected_populations"], list)
    assert set(response["folds"]) == set(PANEL)


def test_cli_simulate_no_perturb_has_no_response_bundle(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "28")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert "response" not in parsed


def test_cli_calibrate_emits_map_and_confidence_intervals(bilan_file: Path) -> None:
    result = _run_cli("calibrate", str(bilan_file), "--bootstrap", "5", "--seed", "1")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert set(parsed["map"]["growth"]) == set(PANEL)
    for pop in PANEL:
        lo, hi = parsed["confidence_intervals"][pop]
        assert lo <= hi
        assert parsed["map"]["growth"][pop] > 0
    assert parsed["n_bootstrap"] == 5
    assert parsed["converged"] is True
    assert isinstance(parsed["mse"], float)


def test_cli_calibrate_writes_params_file(bilan_file: Path, tmp_path: Path) -> None:
    out = tmp_path / "params.json"
    result = _run_cli("calibrate", str(bilan_file), "--bootstrap", "5", "--seed", "1", "--out", str(out))
    assert result.returncode == 0
    assert out.is_file()


def test_cli_simulate_with_calibrated_params_is_consumable(bilan_file: Path, tmp_path: Path) -> None:
    out = tmp_path / "params.json"
    _run_cli("calibrate", str(bilan_file), "--bootstrap", "3", "--seed", "1", "--out", str(out))
    result = _run_cli("simulate", str(bilan_file), "--horizon", "10", "--params", str(out))
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    assert set(parsed["populations"]) == set(PANEL)
    cd8 = parsed["populations"]["CD8"]
    assert len(cd8["times"]) == len(cd8["values"])


def test_cli_simulate_with_ci_emits_confidence_bands(bilan_file: Path, tmp_path: Path) -> None:
    out = tmp_path / "params.json"
    _run_cli("calibrate", str(bilan_file), "--bootstrap", "3", "--seed", "1", "--out", str(out))
    result = _run_cli("simulate", str(bilan_file), "--horizon", "5", "--params", str(out), "--with-ci")
    assert result.returncode == 0
    parsed = json.loads(result.stdout)
    cd8 = parsed["populations"]["CD8"]
    assert "ci_lo" in cd8 and "ci_hi" in cd8
    assert len(cd8["ci_lo"]) == len(cd8["times"])
    assert all(lo <= hi for lo, hi in zip(cd8["ci_lo"], cd8["ci_hi"]))


def test_cli_full_pipeline_bilan_then_calibrate_then_simulate_then_score(bilan_file: Path, tmp_path: Path) -> None:
    """US #25 — drive the whole chain (bilan → calibration → simulation → score) via the CLI."""
    read_result = _run_cli("bilan", str(bilan_file))
    assert read_result.returncode == 0
    assert set(json.loads(read_result.stdout)["populations"]) == set(PANEL)

    params = tmp_path / "params.json"
    cal_result = _run_cli("calibrate", str(bilan_file), "--bootstrap", "3", "--seed", "1",
                          "--out", str(params))
    assert cal_result.returncode == 0
    assert params.is_file()

    sim_result = _run_cli("simulate", str(bilan_file), "--horizon", "10",
                          "--params", str(params), "--perturb", "anti-TNF")
    assert sim_result.returncode == 0
    parsed = json.loads(sim_result.stdout)
    assert set(parsed["perturbed"]["populations"]) == set(PANEL)
    assert "response" in parsed
    score = parsed["response"]["score"]
    assert isinstance(score, float)
    assert parsed["response"]["therapy"] == "anti-TNF"


def test_cli_what_if_compares_scores_across_therapies(bilan_file: Path) -> None:
    """US #27 — explore treatment alternatives by scoring each therapy via the CLI."""
    scores = {}
    for therapy in ("anti-PD1", "anti-TNF", "corticoide"):
        result = _run_cli("simulate", str(bilan_file), "--horizon", "20",
                          "--perturb", therapy)
        assert result.returncode == 0
        parsed = json.loads(result.stdout)
        assert parsed["perturbation"] == therapy
        scores[therapy] = parsed["response"]["score"]
    assert set(scores) == {"anti-PD1", "anti-TNF", "corticoide"}
    assert all(isinstance(v, float) for v in scores.values())
    # a clinician must be able to rank the alternatives — the scores must carry
    # real signal, i.e. the therapies must not all produce the identical number
    assert min(scores.values()) != max(scores.values())


def test_cli_bilan_emits_csv_with_columns_and_rows(bilan_file: Path) -> None:
    result = _run_cli("bilan", str(bilan_file), "--format", "csv")
    assert result.returncode == 0
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "population,value"
    assert len(lines) == len(PANEL) + 1
    assert "CD8,0.15" in lines[1]


def test_cli_simulate_emits_long_format_csv(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "2", "--format", "csv")
    assert result.returncode == 0
    lines = result.stdout.strip().splitlines()
    header, rows = lines[0], lines[1:]
    assert header == "time,scenario,population,value"
    assert rows[0].startswith("0.0,baseline,CD8,")
    assert len(rows) == 200 * len(PANEL)  # 200 grid points over the horizon
    assert all("baseline" in row for row in rows)


def test_cli_simulate_csv_covers_unperturbed_and_perturbed_scenarios(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "2", "--format", "csv",
                      "--perturb", "anti-TNF")
    assert result.returncode == 0
    # only the trajectory table rows carry a scenario (times start at 0.0)
    scenarios = {
        row.split(",")[1]
        for row in result.stdout.splitlines()
        if row.startswith("0.0,")
    }
    assert scenarios == {"unperturbed", "perturbed"}


def test_cli_simulate_csv_perturb_also_emits_comparison_and_score(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "2", "--format", "csv",
                      "--perturb", "anti-PD1")
    assert result.returncode == 0
    assert "population,fold_change,direction" in result.stdout
    header = ("therapy,score,interpretation,reference_score,"
              "is_differential,unexpected_populations")
    assert header in result.stdout


def test_cli_simulate_csv_perturb_emits_perturbation_metadata(bilan_file: Path) -> None:
    """US #21/#22 — the CSV what-if carries the same treatment identity as JSON."""
    result = _run_cli("simulate", str(bilan_file), "--horizon", "2", "--format", "csv",
                      "--perturb", "anti-PD1")
    assert result.returncode == 0
    assert "key,value" in result.stdout
    assert "perturbation,anti-PD1" in result.stdout
    assert "parameter_changes:growth:CD8,+40%" in result.stdout


def test_cli_simulate_csv_no_perturb_has_no_score_table(bilan_file: Path) -> None:
    result = _run_cli("simulate", str(bilan_file), "--horizon", "2", "--format", "csv")
    assert result.returncode == 0
    assert "therapy,score" not in result.stdout


def test_cli_calibrate_emits_csv_with_estimate_and_intervals(bilan_file: Path) -> None:
    result = _run_cli("calibrate", str(bilan_file), "--bootstrap", "3", "--seed", "1",
                      "--format", "csv")
    assert result.returncode == 0
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "population,growth,ci_lo,ci_hi"
    cd8 = [col for col in lines[1].split(",")]
    assert cd8[0] == "CD8"
    assert float(cd8[1]) > 0.0
    # a second metadata table carries the calibration quality signal
    assert "key,value" in result.stdout
    assert "n_bootstrap,3" in result.stdout
    assert "converged," in result.stdout


def test_cli_rejects_unknown_format(bilan_file: Path) -> None:
    result = _run_cli("bilan", str(bilan_file), "--format", "xml")
    assert result.returncode == 2
    assert "format" in result.stderr.lower()