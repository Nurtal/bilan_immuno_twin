"""Unit tests for the simulation module (solver stability, steady-state, horizon)."""

from __future__ import annotations

import numpy as np
import pytest

from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.model import rhs
from bilan_immuno_twin.simulation import (
    default_initial_state,
    simulate,
)

BILAN = {
    "CD8": 0.15,
    "Th1": 0.04,
    "Th2": 0.06,
    "Th17": 0.03,
    "B": 0.12,
    "NK": 0.09,
    "Treg": 0.02,
    "Monocytes": 0.11,
}


def test_simulate_returns_finite_trajectory_over_horizon() -> None:
    params = KineticParameters.defaults()
    x0 = default_initial_state(BILAN)
    trj = simulate(params, x0, horizon=28.0)
    assert trj.times[0] == 0.0
    assert pytest.approx(trj.times[-1], abs=1e-9) == 28.0
    assert trj.values.shape[0] == len(POPULATIONS)
    assert trj.values.shape[1] == trj.times.shape[0]
    assert np.isfinite(trj.values).all()


def test_simulate_keeps_populations_non_negative() -> None:
    params = KineticParameters.defaults()
    x0 = default_initial_state(BILAN)
    trj = simulate(params, x0, horizon=60.0, adaptive=True)
    assert (trj.values >= 0).all()


def test_simulate_radau_handles_stiff_system() -> None:
    # Stiff parameterisation: fast-decaying population (hours) beside slow ones (days)
    params = KineticParameters.defaults()
    params.death = {p: (10.0 if p == "Th17" else v) for p, v in params.death.items()}
    x0 = default_initial_state(BILAN)
    trj = simulate(params, x0, horizon=28.0, method="Radau")
    assert np.isfinite(trj.values).all()
    assert (trj.values >= 0).all()


def test_adaptive_detects_steady_state_and_stops_early() -> None:
    params = KineticParameters.defaults()
    x0 = default_initial_state(BILAN)
    trj = simulate(params, x0, horizon=400.0, adaptive=True)
    assert trj.stopped_early
    assert trj.times[-1] < 400.0
    # at the end of the run the state is a steady state
    final = trj.values[:, -1]
    assert float(np.max(np.abs(rhs(0.0, final, params)))) < 1e-3


def test_adaptive_run_from_steady_state_returns_immediately() -> None:
    params = KineticParameters.defaults()
    x0 = default_initial_state(BILAN)
    # first find the steady state
    trj = simulate(params, x0, horizon=400.0, adaptive=True)
    steady = trj.values[:, -1]
    trj2 = simulate(params, steady, horizon=400.0, adaptive=True)
    assert trj2.times.shape[0] == 1
    assert trj2.stopped_early


def test_non_adaptive_run_does_not_stop() -> None:
    params = KineticParameters.defaults()
    x0 = default_initial_state(BILAN)
    trj = simulate(params, x0, horizon=100.0, adaptive=False)
    assert not trj.stopped_early
    assert pytest.approx(trj.times[-1], abs=1e-9) == 100.0


def test_negative_initial_state_raises() -> None:
    params = KineticParameters.defaults()
    x0 = default_initial_state(BILAN)
    x0[0] = -0.1
    with pytest.raises(ValueError, match="non-negative"):
        simulate(params, x0, adaptive=False)


def test_simulate_is_reproducible() -> None:
    params = KineticParameters.defaults()
    x0 = default_initial_state(BILAN)
    a = simulate(params, x0, horizon=28.0)
    b = simulate(params, x0, horizon=28.0)
    assert np.allclose(a.values, b.values)