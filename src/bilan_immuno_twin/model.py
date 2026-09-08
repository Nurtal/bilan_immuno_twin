"""ODE right-hand side of the graphe de régulations.

The state vector holds the immune populations (variables d'état). At each
timestep the quasi-static cytokine levels (ADR-0005) are computed algebraically
from the populations, and the derivative of each population is:

    dx_i/dt = x_i * ( growth_i * (basal_i + activation_i) * (1 - x_i / capacity_i)
                      - death_i - suppression_i )

where activation_i and suppression_i are weighted sums of Hill functions of the
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


def cytokines_from_state(x: np.ndarray) -> dict[str, float]:
    """Compute quasi-static cytokine levels from the population state."""
    state = {pop: float(x[i]) for i, pop in enumerate(POPULATIONS)}
    return {
        cytokine: sum(coef * state[source] for source, coef in sources.items())
        for cytokine, sources in CYTOKINE_SOURCES.items()
    }


def rhs(t: float, x: np.ndarray, params: KineticParameters) -> np.ndarray:
    """Time derivative of the population state under the given parameters."""
    x = np.asarray(x, dtype=float)
    x = np.clip(x, 0.0, None)

    node_values: dict[str, float] = dict(zip(POPULATIONS, x.tolist()))
    node_values.update(cytokines_from_state(x))

    activation: dict[str, float] = {}
    for target in POPULATIONS:
        total = 0.0
        for source, k in activation_edges_for(target):
            total += params.weights[(source, target, "act")] * hill(node_values[source], k)
        activation[target] = total

    suppression: dict[str, float] = {}
    for target in POPULATIONS:
        total = 0.0
        for source, k in suppression_edges_for(target):
            total += params.weights[(source, target, "sup")] * hill(node_values[source], k)
        suppression[target] = total

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