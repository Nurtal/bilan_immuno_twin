"""T6 — self-consistency validation and regression fixtures."""

import json
from pathlib import Path

from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.simulation import default_initial_state, simulate
from bilan_immuno_twin.validate import (
    DEFAULT_TOL_REL,
    compare_snapshot,
    regression_snapshot,
    self_consistency,
)

FIXTURE = Path(__file__).parent / "fixtures" / "regression.json"
REL_TOL = 1e-4


def _reference() -> dict:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def test_self_consistency_passes_on_synthetic_data() -> None:
    result = self_consistency(seed=0, n_bootstrap=50)
    assert result.passed
    assert result.max_relative_error <= DEFAULT_TOL_REL
    assert result.mean_relative_error < 0.15
    assert result.ci_coverage >= 6


def test_self_consistency_is_deterministic() -> None:
    first = self_consistency(seed=3, n_bootstrap=40)
    second = self_consistency(seed=3, n_bootstrap=40)
    assert first.passed == second.passed
    assert all(
        first.per_population[pop] == second.per_population[pop]
        for pop in POPULATIONS
    )


def test_self_consistency_fails_when_synthetic_noise_is_large() -> None:
    # With noise far above what the calibration can absorb, recovery fails.
    result = self_consistency(seed=0, n_bootstrap=20, noise_sd=0.6, tol_rel=0.05)
    assert not result.passed


def test_self_consistency_detects_planted_truth() -> None:
    """A planted truth far from the prior must recover the planted values."""
    planted = self_consistency(seed=0, n_bootstrap=30)
    for pop in POPULATIONS:
        assert abs(planted.estimated_growth[pop] - planted.true_growth[pop]) <= 0.2


def test_regression_fixture_is_current() -> None:
    snapshot = regression_snapshot()
    drift = compare_snapshot(_reference(), snapshot, rel_tol=REL_TOL)
    assert not drift["simulation"], drift["simulation"]
    assert not drift["calibration"], drift["calibration"]


def test_planted_model_change_is_detected_by_fixture() -> None:
    """A params/model change large enough to matter must trip the fixture."""
    baseline = regression_snapshot()

    # Plant a 10 % model change on a single activation weight, then rerun the
    # exact same simulation protocol as the snapshot.
    changed = KineticParameters.defaults()
    changed.weights[("Th1", "CD8", "act")] *= 1.1
    x0 = default_initial_state(baseline["simulation"]["bilan"])
    trj = simulate(changed, x0, horizon=baseline["simulation"]["horizon"])
    candidate = {
        "simulation": {**baseline["simulation"], "final_populations": {
            pop: float(trj.values[i, -1]) for i, pop in enumerate(POPULATIONS)
        }},
        "calibration": baseline["calibration"],
    }
    drift = compare_snapshot(baseline, candidate, rel_tol=REL_TOL)
    assert drift["simulation"], "planted model change not caught by fixture"


def test_cli_validate_subcommand(tmp_path):
    from bilan_immuno_twin.cli import main

    rc = main(["validate", "--seed", "0", "--bootstrap", "30"])
    assert rc == 0


def test_cli_validate_regression_check_passes(tmp_path):
    from bilan_immuno_twin.cli import main

    rc = main(["validate", "--seed", "0", "--bootstrap", "20",
               "--regression-fixture", str(FIXTURE)])
    assert rc == 0


def test_cli_validate_regression_check_fails_on_drift(tmp_path, capsys):
    from bilan_immuno_twin.cli import main

    fixture = tmp_path / "stale.json"
    fixture.write_text(json.dumps({"calibration": _reference()["calibration"]}), encoding="utf-8")
    rc = main(["validate", "--seed", "0", "--bootstrap", "20",
               "--regression-fixture", str(fixture)])
    assert rc == 3
    out = json.loads(capsys.readouterr().out)
    assert out["regression_check"]["passed"] is False