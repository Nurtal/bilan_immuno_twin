"""Perturbations as modifications of kinetic parameters (ADR-0003).

An immunotherapy is a documented mapping ``perturbation -> {parameter:
relative_delta}``. A relative delta of ``0.4`` means the parameter is
multiplied by ``1.4`` (+40%); ``-0.6`` means ``0.4`` (-60%). Parameters are
addressed by stable identifiers:

- ``growth:<population>``       basal proliferation rate
- ``death:<population>``        basal death rate
- ``weight:<source>-><target>:act|sup``  an edge weight of the graphe
- ``capacity:<population>``     carrying capacity
- ``basal:<population>``        basal activation
"""

from __future__ import annotations

import copy
import re

from bilan_immuno_twin.graph import CYTOKINES, KineticParameters, POPULATIONS

_EDGE_PATTERN = re.compile(r"^weight:(?P<source>\w+)->(?P<target>\w+):(?P<kind>act|sup)$")

# --- Documented perturbation catalogue --------------------------------------
#
# Each entry is interpretable as a change of kinetic parameters:
#   - anti-PD1  : PD-1 blockade relieves CD8+ T cell exhaustion -> CD8+ T
#                 proliferation and IFN-γ-driven activation increase.
#   - anti-TNF  : TNF blockade removes TNF-α-driven monocyte/macrophage
#                 activation -> the TNFa -> Monocytes edge weakens.
#   - corticoïde: corticosteroids depress multiple pro-inflammatory signals
#                 (Th1/Th17/monocyte growth and TNF-α/IL-6 signalling).
PERTURBATION_CATALOGUE: dict[str, dict[str, float]] = {
    "anti-PD1": {
        "growth:CD8": 0.40,
        "death:CD8": -0.10,
        "weight:IFNg->CD8:act": 0.30,
    },
    "anti-TNF": {
        "weight:TNFa->Monocytes:act": -0.60,
    },
    "corticoide": {
        "growth:Th1": -0.35,
        "growth:Th17": -0.40,
        "growth:Monocytes": -0.30,
        "weight:TNFa->Monocytes:act": -0.50,
        "weight:IL6->Th17:act": -0.40,
        "weight:IL6->Monocytes:act": -0.30,
    },
}


class PerturbationError(ValueError):
    """Raised for an unknown therapy or an invalid parameter identifier."""


def known_perturbations() -> list[str]:
    return list(PERTURBATION_CATALOGUE)


def get_perturbation(name: str) -> dict[str, float]:
    """Return the parameter deltas for a named therapy."""
    if name not in PERTURBATION_CATALOGUE:
        known = ", ".join(sorted(PERTURBATION_CATALOGUE))
        raise PerturbationError(f"unknown perturbation '{name}'; known: {known}")
    return dict(PERTURBATION_CATALOGUE[name])


def set_parameter(params: KineticParameters, identifier: str, value: float) -> None:
    """Set a parameter at ``identifier`` to ``value`` (in place)."""
    if identifier.startswith("growth:"):
        pop = identifier[len("growth:"):]
        _require_population(pop, identifier)
        params.growth[pop] = value
        return
    if identifier.startswith("death:"):
        pop = identifier[len("death:"):]
        _require_population(pop, identifier)
        params.death[pop] = value
        return
    if identifier.startswith("capacity:"):
        pop = identifier[len("capacity:"):]
        _require_population(pop, identifier)
        params.capacity[pop] = value
        return
    if identifier.startswith("basal:"):
        pop = identifier[len("basal:"):]
        _require_population(pop, identifier)
        params.basal[pop] = value
        return
    match = _EDGE_PATTERN.match(identifier)
    if match:
        source, target, kind = match.group("source"), match.group("target"), match.group("kind")
        _require_node(source, identifier)
        _require_population(target, identifier)
        params.weights[(source, target, kind)] = value
        return
    raise PerturbationError(f"invalid parameter identifier '{identifier}'")


def apply_perturbation(
    params: KineticParameters, perturbation: dict[str, float]
) -> KineticParameters:
    """Return a copy of ``params`` with every relative delta applied (ADR-0003)."""
    result = copy.deepcopy(params)
    for identifier, delta in perturbation.items():
        if not isinstance(delta, (int, float)):
            raise PerturbationError(
                f"perturbation delta for '{identifier}' must be numeric, got {delta!r}"
            )
        if delta < -1.0:
            raise PerturbationError(
                f"perturbation delta for '{identifier}' is {delta}; "
                "it cannot reduce the parameter below zero"
            )
        current = get_parameter(result, identifier)
        set_parameter(result, identifier, current * (1.0 + delta))
    return result


def get_parameter(params: KineticParameters, identifier: str) -> float:
    """Read the value of a kinetic parameter by its stable identifier."""
    if identifier.startswith("growth:"):
        pop = identifier[len("growth:"):]
        _require_population(pop, identifier)
        return params.growth[pop]
    if identifier.startswith("death:"):
        pop = identifier[len("death:"):]
        _require_population(pop, identifier)
        return params.death[pop]
    if identifier.startswith("capacity:"):
        pop = identifier[len("capacity:"):]
        _require_population(pop, identifier)
        return params.capacity[pop]
    if identifier.startswith("basal:"):
        pop = identifier[len("basal:"):]
        _require_population(pop, identifier)
        return params.basal[pop]
    match = _EDGE_PATTERN.match(identifier)
    if match:
        return params.weights[(match.group("source"), match.group("target"), match.group("kind"))]
    raise PerturbationError(f"invalid parameter identifier '{identifier}'")


def describe(changes: dict[str, float]) -> dict[str, str]:
    """Human-readable interpretation of a set of parameter changes."""
    return {identifier: f"{delta * 100:+.0f}%" for identifier, delta in changes.items()}


def _require_population(name: str, identifier: str) -> None:
    if name not in POPULATIONS:
        raise PerturbationError(
            f"'{identifier}': '{name}' is not a population of the panel "
            f"(known: {', '.join(POPULATIONS)})"
        )


def _require_node(name: str, identifier: str) -> None:
    if name not in POPULATIONS + CYTOKINES:
        raise PerturbationError(
            f"'{identifier}': '{name}' is not a node of the graphe "
            f"(known: {', '.join(POPULATIONS + CYTOKINES)})"
        )