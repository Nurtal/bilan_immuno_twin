"""ODE right-hand side of the graphe de régulations.

The state vector holds the immune populations (variables d'état). At each
timestep the quasi-static cytokine levels (ADR-0005) are computed algebraically
from the populations, and the derivative of each population is:

    dx_i/dt = x_i * ( growth_i * (basal_i + activation_i) * damp_i - death_i )

with the logistic saturation `(1 - x_i / capacity_i)` folded in, where
activation_i and suppression_i are weighted sums of Hill functions of the
regulatory node values (populations or cytokines). The x_i prefactor gives
zero-stability: an absent population stays absent (no spontaneous generation),
keeping the solution non-negative for non-negative initial states.
"""

from __future__ import annotations

import numpy as np

from bilan_immuno_twin.graph import (
    CYTOKINE_SOURCES,
    KineticParameters,
    POPULATIONS,
    activation_edges_for,
    suppression_edges_for,
)
from bilan_immuno_twin.hill import hill

STEADY_STATE_MAX_ITER = 500
STEADY_STATE_TOL = 1e-12


def cytokines_from_state(x: np.ndarray) -> dict[str, float]:
    """Compute quasi-static cytokine levels from the population state."""
    state = {pop: float(x[i]) for i, pop in enumerate(POPULATIONS)}
    return {
        cytokine: sum(coef * state[source] for source, coef in sources.items())
        for cytokine, sources in CYTOKINE_SOURCES.items()
    }


def node_values_from_state(x: np.ndarray) -> dict[str, float]:
    """All graph node values (populations + quasi-static cytokines)."""
    node_values: dict[str, float] = dict(zip(POPULATIONS, x.tolist()))
    node_values.update(cytokines_from_state(x))
    return node_values


def activation_and_suppression(
    node_values: dict[str, float], params: KineticParameters
) -> tuple[dict[str, float], dict[str, float]]:
    """Weighted-Hill activation and suppression input per population."""
    activation: dict[str, float] = {}
    for target in POPULATIONS:
        activation[target] = sum(
            params.weights[(source, target, "act")] * hill(node_values[source], k)
            for source, k in activation_edges_for(target)
        )
    suppression: dict[str, float] = {}
    for target in POPULATIONS:
        suppression[target] = sum(
            params.weights[(source, target, "sup")] * hill(node_values[source], k)
            for source, k in suppression_edges_for(target)
        )
    return activation, suppression


def rhs(t: float, x: np.ndarray, params: KineticParameters) -> np.ndarray:
    """Time derivative of the population state under the given parameters."""
    x = np.clip(np.asarray(x, dtype=float), 0.0, None)
    activation, suppression = activation_and_suppression(
        node_values_from_state(x), params
    )

    deriv = np.empty(len(POPULATIONS), dtype=float)
    for i, pop in enumerate(POPULATIONS):
        x_pop = x[i]
        logistic = (1.0 - x_pop / params.capacity[pop])
        if logistic < 0.0:
            logistic = 0.0
        # suppression is saturating and multiplicative: it damps growth but
        # cannot by itself drive a population extinct (ADR-0001 spirit)
        damp = 1.0 / (1.0 + suppression[pop])
        growth = params.growth[pop] * (params.basal[pop] + activation[pop]) * damp
        deriv[i] = x_pop * (growth * logistic - params.death[pop])
    return deriv


def steady_state(
    params: KineticParameters,
    x0: np.ndarray,
    max_iter: int = STEADY_STATE_MAX_ITER,
    tol: float = STEADY_STATE_TOL,
    x_tol: float = 1e-12,
) -> np.ndarray:
    """Model steady state reached from ``x0``, by fixed-point iteration.

    At a nonzero equilibrium ``dx_i/dt = 0`` gives a closed form for each
    population in terms of the others:

        x_i = capacity_i * (1 - death_i / (growth_i * (basal_i + act_i) * damp_i))

    so the steady state is found by iterating this map to a fixed point. Cheap
    and smooth in the parameters, which is what the calibrator needs.
    """
    x = np.clip(np.asarray(x0, dtype=float), 0.0, None)
    for _ in range(max_iter):
        activation, suppression = activation_and_suppression(
            node_values_from_state(x), params
        )
        x_new = np.empty(len(POPULATIONS), dtype=float)
        for i, pop in enumerate(POPULATIONS):
            damp = 1.0 / (1.0 + suppression[pop])
            growth_total = params.growth[pop] * (params.basal[pop] + activation[pop]) * damp
            capacity = params.capacity[pop]
            if growth_total > params.death[pop]:
                x_new[i] = capacity * (1.0 - params.death[pop] / growth_total)
            else:
                x_new[i] = 0.0
            x_new[i] = min(x_new[i], capacity)
        if np.max(np.abs(x_new - x)) < x_tol:
            return np.clip(x_new, 0.0, None)
        x = x_new
    return np.clip(x, 0.0, None)