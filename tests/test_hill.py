"""Unit tests for the fonctions de Hill (pure-math seam, ADR-0001)."""

from __future__ import annotations

import math

import pytest

from bilan_immuno_twin.hill import hill


def test_hill_zero_is_zero() -> None:
    assert hill(0.0, k=0.1) == 0.0


def test_hill_at_threshold_is_half() -> None:
    assert math.isclose(hill(0.2, k=0.2), 0.5)


def test_hill_saturates_to_one() -> None:
    assert math.isclose(hill(1000.0, k=0.2), 1.0, rel_tol=1e-6)


def test_hill_monotonically_increasing() -> None:
    x = [0.01, 0.05, 0.2, 0.8, 3.0, 10.0]
    values = [hill(v, k=0.2, n=2) for v in x]
    assert all(a < b for a, b in zip(values, values[1:]))


def test_hill_supra_threshold_response_is_sigmoid() -> None:
    # below K -> small response; well above -> near saturation
    assert hill(0.01, k=0.2, n=2) < 0.01
    assert hill(0.2, k=0.2, n=2) == 0.5
    assert hill(2.0, k=0.2, n=2) > 0.99


@pytest.mark.parametrize("n", [1, 2, 4])
def test_hill_higher_cooperativity_sharpens_threshold(n: float) -> None:
    # at x slightly below K the response is smaller for steeper n
    assert hill(0.5, k=1.0, n=4) < hill(0.5, k=1.0, n=1)


def test_hill_negative_x_returns_zero() -> None:
    assert hill(-5.0, k=0.2) == 0.0


def test_hill_invalid_k_raises() -> None:
    with pytest.raises(ValueError):
        hill(0.5, k=0.0)


def test_hill_invalid_n_raises() -> None:
    with pytest.raises(ValueError):
        hill(0.5, k=0.2, n=0.0)