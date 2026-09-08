"""Unit tests for the ODE right-hand side (pure-math seam, T2)."""

from __future__ import annotations

import math

import numpy as np

from bilan_immuno_twin.graph import CYTOKINE_SOURCES, KineticParameters, POPULATIONS
from bilan_immuno_twin.model import cytokines_from_state, rhs

# A mid-range immune state (fractions of total immune cells).
BILAN = np.array([0.15, 0.04, 0.06, 0.03, 0.12, 0.09, 0.02, 0.11], dtype=float)


def test_rhs_returns_same_dimension() -> None:
    params = KineticParameters.defaults()
    deriv = rhs(0.0, BILAN, params)
    assert deriv.shape == BILAN.shape
    assert not np.isnan(deriv).any()


def test_zero_stability_absent_population_stays_absent() -> None:
    # if a population is 0, its derivative must be 0 (no spontaneous generation)
    params = KineticParameters.defaults()
    x = BILAN.copy()
    x[0] = 0.0  # no CD8
    deriv = rhs(0.0, x, params)
    assert deriv[0] == 0.0
    # also at the zero state everything stays zero
    assert (rhs(0.0, np.zeros(len(POPULATIONS)), params) == 0.0).all()


def test_rhs_negative_state_is_clipped() -> None:
    params = KineticParameters.defaults()
    x = BILAN.copy()
    x[2] = -0.5
    deriv = rhs(0.0, x, params)
    assert not np.isnan(deriv).any()


def test_beyond_capacity_population_is_pulled_back() -> None:
    # a population far above its carrying capacity should decay (negative derivative)
    params = KineticParameters.defaults()
    x = BILAN.copy()
    x[0] = 3.0  # CD8 far above capacity 1.0
    deriv = rhs(0.0, x, params)
    assert deriv[0] < 0.0


def test_stronger_activation_increases_population_growth() -> None:
    # raising an activation weight increases the corresponding derivative when
    # the activator is present
    params = KineticParameters.defaults()
    x = BILAN.copy()
    base = rhs(0.0, x, params)
    params.weights[("Th1", "CD8", "act")] = 5.0
    boosted = rhs(0.0, x, params)
    assert boosted[0] > base[0]


def test_stronger_suppression_decreases_population_growth() -> None:
    params = KineticParameters.defaults()
    x = BILAN.copy()
    base = rhs(0.0, x, params)
    params.weights[("Treg", "CD8", "sup")] = 5.0
    suppressed = rhs(0.0, x, params)
    assert suppressed[0] < base[0]


def test_zero_growth_with_no_activation_and_no_basal_is_decay() -> None:
    # with basal activation removed and no activators present, only death remains
    params = KineticParameters.defaults()
    params.basal = {p: 0.0 for p in POPULATIONS}
    params.weights = {k: 0.0 for k in params.weights}
    x = BILAN.copy()
    x[0] = 0.5  # isolated population with no input
    deriv = rhs(0.0, x, params)
    assert math.isclose(deriv[0], -params.death["CD8"] * 0.5)


def test_cytokines_are_algebraic_functions_of_state() -> None:
    values = cytokines_from_state(BILAN)
    assert set(values) == set(CYTOKINE_SOURCES)
    for cytokine, sources in CYTOKINE_SOURCES.items():
        expected = sum(coef * float(BILAN[POPULATIONS.index(pop)])
                       for pop, coef in sources.items())
        assert math.isclose(values[cytokine], expected)


def test_cytokines_of_zero_state_are_zero() -> None:
    assert all(v == 0.0 for v in cytokines_from_state(np.zeros(len(POPULATIONS))).values())


def test_rhs_time_independent() -> None:
    params = KineticParameters.defaults()
    assert np.allclose(rhs(0.0, BILAN, params), rhs(42.0, BILAN, params))