"""Simulation of the graphe de régulations ODE system.

Integrates the system with `scipy.integrate.solve_ivp` using the Radau
implicit method to handle stiffness (raideur). Supports a fixed horizon and an
adaptive mode that stops early once the system reaches a steady state.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.model import rhs

STEADY_STATE_ABS_TOL = 1e-4
STEADY_STATE_METHOD = "Radau"


@dataclass
class Trajectory:
    times: np.ndarray  # shape (n_times,)
    values: np.ndarray  # shape (n_populations, n_times)
    stopped_early: bool = False
    message: str = ""


def default_initial_state(bilan: dict) -> np.ndarray:
    """Ordered initial populations from a parsed bilan."""
    return np.array([bilan[pop] for pop in POPULATIONS], dtype=float)


def simulate(
    params: KineticParameters,
    x0: np.ndarray,
    horizon: float = 28.0,
    n_times: int = 200,
    adaptive: bool = False,
    method: str = STEADY_STATE_METHOD,
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> Trajectory:
    """Simulate the immune population dynamics over [0, horizon] days.

    Integrates with an implicit solver (Radau by default) to handle stiffness.
    With ``adaptive=True``, the run terminates as soon as every population's
    derivative magnitude drops below ``STEADY_STATE_ABS_TOL``, and the returned
    trajectory ends at the steady-state time.
    """
    x0 = np.asarray(x0, dtype=float)
    if np.any(x0 < 0):
        raise ValueError("initial state must be non-negative")

    if adaptive and _max_derivative(x0, params) < STEADY_STATE_ABS_TOL:
        return Trajectory(
            times=np.array([0.0]),
            values=x0.reshape(-1, 1),
            stopped_early=True,
            message="initial state is already a steady state",
        )

    t_span = (0.0, float(horizon))
    events = None
    if adaptive:
        events = [_make_steady_state_event(params, STEADY_STATE_ABS_TOL)]

    sol = solve_ivp(
        lambda t, x: rhs(t, x, params),
        t_span,
        y0=x0,
        method=method,
        t_eval=np.linspace(t_span[0], t_span[1], max(n_times, 2)),
        rtol=rtol,
        atol=atol,
        events=events,
    )

    stopped_early = bool(sol.t_events is not None and any(len(e) > 0 for e in sol.t_events))
    return Trajectory(
        times=sol.t,
        # clip tiny solver negatives so states stay physically non-negative
        values=np.clip(sol.y, 0.0, None),
        stopped_early=stopped_early,
        message=str(sol.message),
    )


def _max_derivative(x: np.ndarray, params: KineticParameters) -> float:
    return float(np.max(np.abs(rhs(0.0, x, params))))


def _make_steady_state_event(params: KineticParameters, tol: float):
    """Terminal event fired (direction -1) when the system settles to steady state."""
    def event(t: float, x: np.ndarray) -> float:
        return float(np.max(np.abs(rhs(t, x, params)))) - tol

    event.terminal = True
    event.direction = -1
    return event