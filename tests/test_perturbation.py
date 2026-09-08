"""Unit tests for perturbation mapping (ADR-0003, T3 pure seam)."""

from __future__ import annotations

import pytest

from bilan_immuno_twin.graph import KineticParameters
from bilan_immuno_twin.perturbation import (
    PerturbationError,
    apply_perturbation,
    get_perturbation,
    get_parameter,
    known_perturbations,
    set_parameter,
)


def test_catalogue_contains_opening_therapies() -> None:
    assert {"anti-PD1", "anti-TNF", "corticoide"} <= set(known_perturbations())


def test_anti_pd1_raises_cd8_activation() -> None:
    changes = get_perturbation("anti-PD1")
    assert changes["growth:CD8"] > 0
    assert changes["weight:IFNg->CD8:act"] > 0


def test_anti_tnf_reduces_tnf_effect() -> None:
    changes = get_perturbation("anti-TNF")
    assert changes["weight:TNFa->Monocytes:act"] < 0


def test_corticoide_depresses_proinflammatory_signals() -> None:
    changes = get_perturbation("corticoide")
    assert changes["growth:Th1"] < 0
    assert changes["growth:Th17"] < 0
    assert changes["growth:Monocytes"] < 0


def test_every_catalogue_identifier_names_a_live_parameter() -> None:
    params = KineticParameters.defaults()
    for therapy in known_perturbations():
        for identifier in get_perturbation(therapy):
            get_parameter(params, identifier)  # must not raise


def test_unknown_therapy_raises_clear_error() -> None:
    with pytest.raises(PerturbationError, match="unknown perturbation"):
        get_perturbation("vax-covid")


def test_accented_corticoïde_spelling_is_accepted() -> None:
    accented = get_perturbation("corticoïde")
    canonical = get_perturbation("corticoide")
    assert accented == canonical


def test_apply_relative_delta_multiplicative() -> None:
    params = KineticParameters.defaults()
    base = params.growth["CD8"]
    out = apply_perturbation(params, {"growth:CD8": 0.40})
    assert out.growth["CD8"] == pytest.approx(base * 1.4)
    # original is unchanged (no mutation)
    assert params.growth["CD8"] == base


def test_apply_negative_delta_reduces_parameter() -> None:
    params = KineticParameters.defaults()
    base = params.weights[("TNFa", "Monocytes", "act")]
    out = apply_perturbation(params, {"weight:TNFa->Monocytes:act": -0.60})
    assert out.weights[("TNFa", "Monocytes", "act")] == pytest.approx(base * 0.4)


def test_delta_below_minus_one_rejected() -> None:
    params = KineticParameters.defaults()
    with pytest.raises(PerturbationError, match="below zero"):
        apply_perturbation(params, {"growth:CD8": -1.5})


def test_invalid_identifier_rejected() -> None:
    params = KineticParameters.defaults()
    with pytest.raises(PerturbationError, match="not a population"):
        apply_perturbation(params, {"growth:Nonsense": 0.1})
    with pytest.raises(PerturbationError, match="invalid parameter identifier"):
        apply_perturbation(params, {"whatever:CD8": 0.1})


def test_set_parameter_supports_and_returns_value() -> None:
    params = KineticParameters.defaults()
    set_parameter(params, "growth:CD8", 1.2)
    assert get_parameter(params, "growth:CD8") == 1.2


def test_catalogue_is_interpretable_as_documented_changes() -> None:
    # every therapy entry must be expressible in the documented vocabulary
    from bilan_immuno_twin.perturbation import describe

    for therapy in known_perturbations():
        description = describe(get_perturbation(therapy))
        for identifier, percent in description.items():
            assert identifier.startswith(("growth:", "death:", "weight:", "capacity:", "basal:"))
            assert percent.endswith("%")