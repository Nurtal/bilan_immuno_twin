"""Unit tests for MAP calibration + bootstrap (T5 pure seam, ADR-0002)."""

from __future__ import annotations

import numpy as np
import pytest

from bilan_immuno_twin.calibration import (
    GROWTH_BOUNDS,
    calibrate,
    parameters_from_growth,
)
from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.model import steady_state

TRUTH_NEAR_PRIOR = dict(zip(POPULATIONS, [0.32, 0.36, 0.33, 0.38, 0.34, 0.31, 0.37, 0.35]))
TRUTH_SPREAD = dict(zip(POPULATIONS, [0.28, 0.42, 0.30, 0.38, 0.40, 0.29, 0.26, 0.44]))

BILAN = np.array([0.15, 0.04, 0.06, 0.03, 0.12, 0.09, 0.02, 0.11])


def _synthetic_bilan(growth: dict[str, float],
                     rng: np.random.Generator | None = None) -> dict[str, float]:
    params = KineticParameters.defaults()
    params.growth = dict(growth)
    ss = steady_state(params, BILAN)
    observed = dict(zip(POPULATIONS, ss.tolist()))
    if rng is not None:
        observed = dict(zip(POPULATIONS, (ss * np.exp(rng.normal(0.0, 0.05, len(POPULATIONS)))).tolist()))
    return observed


def test_map_recovers_known_parameters_on_synthetic_data() -> None:
    observed = _synthetic_bilan(TRUTH_NEAR_PRIOR)
    result = calibrate(observed, n_bootstrap=20, seed=0)
    assert result.converged
    for pop in POPULATIONS:
        assert abs(result.growth_estimates[pop] - TRUTH_NEAR_PRIOR[pop]) / TRUTH_NEAR_PRIOR[pop] < 0.20


def test_bootstrap_intervals_contain_truth() -> None:
    observed = _synthetic_bilan(TRUTH_SPREAD)
    result = calibrate(observed, n_bootstrap=40, seed=2)
    covered = sum(
        result.confidence_intervals[pop][0] <= TRUTH_SPREAD[pop] <= result.confidence_intervals[pop][1]
        for pop in POPULATIONS
    )
    assert covered >= 6


def test_calibration_reproduces_observed_steady_state() -> None:
    observed = _synthetic_bilan(TRUTH_SPREAD)
    result = calibrate(observed, n_bootstrap=10, seed=3)
    params = parameters_from_growth(result.growth_estimates)
    y = np.array([observed[p] for p in POPULATIONS])
    ss = steady_state(params, y)
    # the calibrated model holds the observed state near a steady state
    assert float(np.max(np.abs(np.log(ss) - np.log(y)))) < 0.35


def test_calibration_is_deterministic_with_seed() -> None:
    observed = _synthetic_bilan(TRUTH_NEAR_PRIOR)
    a = calibrate(observed, n_bootstrap=10, seed=42)
    b = calibrate(observed, n_bootstrap=10, seed=42)
    assert a.growth_estimates == b.growth_estimates
    assert np.array_equal(a.bootstrap_values, b.bootstrap_values)


def test_calibrated_parameters_are_consumable_by_simulation() -> None:
    from bilan_immuno_twin.simulation import simulate

    observed = _synthetic_bilan(TRUTH_NEAR_PRIOR)
    result = calibrate(observed, n_bootstrap=5, seed=5)
    params = parameters_from_growth(result.growth_estimates)
    trj = simulate(params, BILAN, horizon=28.0)
    assert np.isfinite(trj.values).all()
    assert (trj.values >= 0).all()


def test_estimates_stay_within_support_bounds() -> None:
    observed = _synthetic_bilan(TRUTH_SPREAD)
    result = calibrate(observed, n_bootstrap=10, seed=6)
    lo, hi = GROWTH_BOUNDS
    for pop in POPULATIONS:
        assert lo <= result.growth_estimates[pop] <= hi
        assert lo <= result.confidence_intervals[pop][0] <= result.confidence_intervals[pop][1] <= hi


def test_growth_rates_changed_vs_defaults() -> None:
    observed = _synthetic_bilan(TRUTH_SPREAD)
    result = calibrate(observed, n_bootstrap=10, seed=7)
    params = parameters_from_growth(result.growth_estimates)
    defaults = KineticParameters.defaults()
    assert params.growth != defaults.growth


def test_negative_populations_rejected() -> None:
    bad = dict(_synthetic_bilan(TRUTH_NEAR_PRIOR))
    bad["CD8"] = -0.1
    with pytest.raises(ValueError, match="non-negative"):
        calibrate(bad, n_bootstrap=1)