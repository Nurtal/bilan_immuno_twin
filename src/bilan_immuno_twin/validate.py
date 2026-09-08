"""Self-consistency validation and regression fixtures (T6).

The framework is validated on synthetic data before any clinical data (out of
scope for the MVP). The self-consistency protocol is:

    fixed parameters -> simulate -> add noise -> recalibrate -> recover

A run passes when every calibrated growth rate is recovered to within a
relative tolerance. Regression snapshots pin deterministic simulation and
calibration outputs so that model changes that silently alter results are
caught.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from bilan_immuno_twin.calibration import calibrate
from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.model import steady_state
from bilan_immuno_twin.response import REFERENCE_BILAN
from bilan_immuno_twin.simulation import default_initial_state, simulate

DEFAULT_NOISE_SD = 0.05     # synthetic measurement noise (log scale) - careful
DEFAULT_TOL_REL = 0.25      # recovery within ±25 % relative
REGRESSION_HORIZON = 28.0

# Truths chosen near the prior so recovery is identifiable.
TRUTH_GROWTH = dict(zip(POPULATIONS, [0.32, 0.36, 0.33, 0.38, 0.34, 0.31, 0.37, 0.35]))


@dataclass
class SelfConsistencyResult:
    passed: bool
    max_relative_error: float
    mean_relative_error: float
    per_population: dict[str, float]  # population -> relative error
    true_growth: dict[str, float]
    estimated_growth: dict[str, float]
    ci_coverage: int
    confidence_intervals: dict[str, tuple[float, float]] = field(default_factory=dict)


def self_consistency(
    seed: int = 0,
    noise_sd: float = DEFAULT_NOISE_SD,
    n_bootstrap: int = 50,
    tol_rel: float = DEFAULT_TOL_REL,
) -> SelfConsistencyResult:
    """Fixed params -> simulate -> noise -> recalibrate -> recover (within tol)."""
    rng = np.random.default_rng(seed)
    params = KineticParameters.defaults()
    params.growth = dict(TRUTH_GROWTH)

    x0 = default_initial_state(REFERENCE_BILAN)
    ss = steady_state(params, x0)
    observed = dict(zip(
        POPULATIONS,
        np.clip(ss * np.exp(rng.normal(0.0, noise_sd, len(POPULATIONS))), 1e-9, None).tolist(),
    ))

    result = calibrate(observed, n_bootstrap=n_bootstrap, seed=seed)

    errors = {
        p: abs(result.growth_estimates[p] - TRUTH_GROWTH[p]) / TRUTH_GROWTH[p]
        for p in POPULATIONS
    }
    max_error = float(np.max(list(errors.values())))
    mean_error = float(np.mean(list(errors.values())))
    coverage = sum(
        result.confidence_intervals[p][0] <= TRUTH_GROWTH[p] <= result.confidence_intervals[p][1]
        for p in POPULATIONS
    )

    return SelfConsistencyResult(
        passed=max_error <= tol_rel,
        max_relative_error=max_error,
        mean_relative_error=mean_error,
        per_population=errors,
        true_growth=dict(TRUTH_GROWTH),
        estimated_growth=result.growth_estimates,
        ci_coverage=coverage,
        confidence_intervals=result.confidence_intervals,
    )


def regression_snapshot() -> dict:
    """Deterministic snapshot of simulation + calibration outputs for a fixed input.

    Used to pin outputs across versions (regression fixtures).
    """
    params = KineticParameters.defaults()
    x0 = default_initial_state(REFERENCE_BILAN)
    trj = simulate(params, x0, horizon=REGRESSION_HORIZON)
    simulation = {
        pop: float(trj.values[i, -1]) for i, pop in enumerate(POPULATIONS)
    }

    result = calibrate(REFERENCE_BILAN, n_bootstrap=20, seed=0)
    calibration = {
        pop: {
            "map": result.growth_estimates[pop],
            "ci_lo": result.confidence_intervals[pop][0],
            "ci_hi": result.confidence_intervals[pop][1],
        }
        for pop in POPULATIONS
    }

    return {
        "simulation": {
            "bilan": dict(REFERENCE_BILAN),
            "horizon": REGRESSION_HORIZON,
            "final_populations": simulation,
        },
        "calibration": {
            "bilan": dict(REFERENCE_BILAN),
            "n_bootstrap": 20,
            "seed": 0,
            "growth": calibration,
        },
    }


def compare_snapshot(reference: dict, current: dict, rel_tol: float = 1e-4) -> dict:
    """Compare two snapshots and report which fields drift beyond tolerance.

    A section or population missing from the reference fixture counts as drift
    (the fixture is stale relative to the current model).
    """
    drift: dict[str, dict[str, float]] = {}
    for section in ("simulation", "calibration"):
        drift[section] = {}
        ref_pops = (reference.get(section) or {}).get(
            "final_populations" if section == "simulation" else "growth", {}
        )
        cur_pops = (current.get(section) or {}).get(
            "final_populations" if section == "simulation" else "growth", {}
        )
        for pop in POPULATIONS:
            if pop not in ref_pops:
                drift[section][pop] = float("inf")
                continue
            ref_value = ref_pops[pop] if isinstance(ref_pops[pop], (int, float)) else ref_pops[pop]["map"]
            cur_value = cur_pops[pop] if isinstance(cur_pops[pop], (int, float)) else cur_pops[pop]["map"]
            relative = abs(cur_value - ref_value) / max(abs(ref_value), 1e-9)
            if relative > rel_tol:
                drift[section][pop] = float(relative)
    return drift