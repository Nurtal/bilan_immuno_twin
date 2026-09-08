"""Calibration bayésienne: MAP + bootstrap (ADR-0002).

A bilan is a single cytometric snapshot, so rates are not directly observable.
The calibration assumes the observed populations are a steady state of the
graphe de régulations. Writing the equilibrium condition ``dx_i/dt = 0`` at the
observed state ``y`` for each nonzero population gives a closed form for the
per-population basal growth rate:

    g*_i = death_i / ( (basal_i + act_i(y)) * damp_i(y) * (1 - y_i / M_i) )

where ``act_i(y)`` and ``damp_i(y)`` depend only on the observed state. The MAP
estimate shrinks ``g*`` toward the literature prior mean in log space:

    log g_map = (log g* / sigma_obs^2 + log mu / sigma_prior^2)
                / (1 / sigma_obs^2 + 1 / sigma_prior^2)

Bootstrap resampling (multiplicative log-normal noise on the observed bilan)
propagates measurement uncertainty into confidence intervals.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.model import activation_and_suppression, node_values_from_state

# Literature-derived priors (log space).
GROWTH_PRIOR_LOG_MEAN = np.log(0.35)     # default growth rate (per day)
GROWTH_PRIOR_LOG_SD = 0.5                # moderate spread across patients
MEASUREMENT_LOG_SD = 0.15                # ~15 % CV on cytometric fractions
GROWTH_BOUNDS = (0.02, 2.0)
_EPS = 1e-9

_LOG_SIGMA_OBS_2 = MEASUREMENT_LOG_SD ** 2
_LOG_SIGMA_PRIOR_2 = GROWTH_PRIOR_LOG_SD ** 2


@dataclass
class CalibrationResult:
    growth_estimates: dict[str, float]
    confidence_intervals: dict[str, tuple[float, float]]
    bootstrap_values: np.ndarray  # shape (n_bootstrap, n_populations)
    converged: bool
    mse: float
    message: str = ""
    infeasible: list[str] = field(default_factory=list)


def _locally_exact_growth(y: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """Growth rates that exactly hold the observed state at equilibrium.

    Returns (g*, infeasible) where g*[i] = death_i / c_i with
    c_i = (basal_i + act_i(y)) * damp_i(y) * (1 - y_i / M_i). Populations for
    which the observed level exceeds their carrying capacity (c_i <= 0) or is
    zero (any growth rate leaves the derivative zero) are infeasible and
    receive the prior mean.
    """
    params = KineticParameters.defaults()
    activation, suppression = activation_and_suppression(node_values_from_state(y), params)

    g_star = np.empty(len(POPULATIONS), dtype=float)
    infeasible: list[str] = []
    for i, pop in enumerate(POPULATIONS):
        y_i = float(y[i])
        if y_i <= _EPS:
            g_star[i] = float(np.exp(GROWTH_PRIOR_LOG_MEAN))
            continue
        damp = 1.0 / (1.0 + suppression[pop])
        logistic = max(0.0, 1.0 - y_i / params.capacity[pop])
        c_i = (params.basal[pop] + activation[pop]) * damp * logistic
        if c_i <= _EPS:
            infeasible.append(pop)
            g_star[i] = float(np.exp(GROWTH_PRIOR_LOG_MEAN))
        else:
            g_star[i] = params.death[pop] / c_i
    return g_star, infeasible


def _posterior_shrink(g_star: np.ndarray, y: np.ndarray) -> np.ndarray:
    """MAP growth from closed-form MLE plus literature prior (log space).

    Data weight is discounted for populations sitting near their carrying
    capacity: there the observed level barely constrains growth, so the prior
    should dominate rather than letting small measurement noise produce extreme
    rate estimates.
    """
    log_gstar = np.log(np.clip(g_star, _EPS, None))
    params = KineticParameters.defaults()
    w_prior = 1.0 / _LOG_SIGMA_PRIOR_2
    log_map = np.empty(len(POPULATIONS), dtype=float)
    for i, _pop in enumerate(POPULATIONS):
        capacity = params.capacity[POPULATIONS[i]]
        proximity = float(y[i]) / max(capacity - float(y[i]), _EPS)
        discount = 1.0 / (1.0 + proximity ** 2)
        w_data = (1.0 / _LOG_SIGMA_OBS_2) * discount
        log_map[i] = (w_data * log_gstar[i] + w_prior * GROWTH_PRIOR_LOG_MEAN) / (w_data + w_prior)
    lo = np.log(GROWTH_BOUNDS[0])
    hi = np.log(GROWTH_BOUNDS[1])
    return np.clip(log_map, lo, hi)


def _resample_observed(observed: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Multiplicative log-normal resample of the observed bilan (bootstrap)."""
    noise = rng.normal(0.0, MEASUREMENT_LOG_SD, size=observed.shape)
    return np.clip(observed * np.exp(noise), _EPS, None)


def _fit_growth(observed: np.ndarray) -> tuple[np.ndarray, list[str]]:
    y = np.clip(np.asarray(observed, dtype=float), _EPS, None)
    g_star, infeasible = _locally_exact_growth(y)
    return _posterior_shrink(g_star, y), infeasible


def calibrate(
    observed: dict[str, float],
    n_bootstrap: int = 100,
    seed: int | None = None,
) -> CalibrationResult:
    """MAP-estimate patient-specific growth rates from a bilan + bootstrap CIs."""
    y = np.array([observed[p] for p in POPULATIONS], dtype=float)
    if np.any(y < 0):
        raise ValueError("observed populations must be non-negative")

    log_map, infeasible = _fit_growth(y)
    growth_estimates = {pop: float(np.exp(log_map[i])) for i, pop in enumerate(POPULATIONS)}

    rng = np.random.default_rng(seed)
    bootstrap_values = np.empty((n_bootstrap, len(POPULATIONS)), dtype=float)
    for b in range(n_bootstrap):
        y_b = _resample_observed(y, rng)
        log_b, _ = _fit_growth(y_b)
        bootstrap_values[b] = np.exp(log_b)

    lo = np.percentile(bootstrap_values, 2.5, axis=0)
    hi = np.percentile(bootstrap_values, 97.5, axis=0)
    confidence_intervals = {
        pop: (float(lo[i]), float(hi[i])) for i, pop in enumerate(POPULATIONS)
    }

    g_star, _ = _locally_exact_growth(y)
    mse = float(np.mean((np.exp(log_map) - g_star) ** 2))

    return CalibrationResult(
        growth_estimates=growth_estimates,
        confidence_intervals=confidence_intervals,
        bootstrap_values=bootstrap_values,
        converged=not infeasible,
        mse=mse,
        infeasible=infeasible,
    )


def parameters_from_growth(growth: dict[str, float]) -> KineticParameters:
    """Build simulation parameters from calibrated growth rates."""
    params = KineticParameters.defaults()
    params.growth = dict(growth)
    return params