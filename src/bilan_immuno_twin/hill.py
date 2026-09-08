"""Fonctions de Hill for regulatory interactions (ADR-0001).

`hill(x) = x^n / (K^n + x^n)` captures biological thresholds and saturation:
a signal must reach a level K to activate a response, and the effect saturates.
"""

from __future__ import annotations


def hill(x: float, k: float, n: float = 2.0) -> float:
    """Hill function with threshold K and cooperativity n.

    Properties (tested at the pure-math seam):
      - hill(0) = 0, hill(x) -> 1 as x -> infinity (saturation)
      - hill(K) = 0.5 (threshold)
      - strictly increasing in x for K > 0, n > 0
    """
    if k <= 0:
        raise ValueError(f"Hill threshold K must be positive, got {k}")
    if n <= 0:
        raise ValueError(f"Hill cooperativity n must be positive, got {n}")
    if x < 0:
        return 0.0
    x_n = x ** n
    return x_n / (k ** n + x_n)