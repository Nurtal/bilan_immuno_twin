"""Unit tests for the score de réponse derivation and differential flagging."""

from __future__ import annotations

from bilan_immuno_twin.response import (
    TARGET_DIRECTIONS,
    compute_score,
    evaluate_response,
    fold_changes,
    score_interpretation,
    unexpected_populations,
)


def test_score_is_zero_when_no_change() -> None:
    folds = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}
    assert compute_score(folds, "anti-PD1") == 0.0


def test_score_increases_with_effector_expansion_for_anti_pd1() -> None:
    folds = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}
    folds["CD8"] = 4.0  # doubling twice
    folds["NK"] = 2.0
    score = compute_score(folds, "anti-PD1")
    assert score > 0.0
    # favorable qualitative band, but still continuous
    assert score_interpretation(score) == "favorable"


def test_score_decreases_for_opposite_response() -> None:
    folds = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}
    folds["CD8"] = 0.25
    score = compute_score(folds, "anti-PD1")
    assert score < 0.0
    assert score_interpretation(score) == "unfavorable"


def test_score_is_continuous_not_binary() -> None:
    # two different response magnitudes give two different scores
    mild = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}
    mild["CD8"] = 1.5
    strong = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}
    strong["CD8"] = 3.0
    assert compute_score(mild, "anti-PD1") < compute_score(strong, "anti-PD1")
    assert compute_score(mild, "anti-PD1") != compute_score(strong, "anti-PD1")


def test_corticoide_rules_are_reversed() -> None:
    # for corticoïdes a pro-inflammatory *decrease* is favorable
    folds = {pop: 1.0 for pop in TARGET_DIRECTIONS["corticoide"]}
    folds["Th17"] = 0.25
    folds["Th1"] = 0.5
    assert compute_score(folds, "corticoide") > 0.0


def test_fold_changes_from_final_values() -> None:
    unperturbed = {"CD8": 1.0, "Th1": 1.0}
    perturbed = {"CD8": 3.0, "Th1": 0.5}
    # tail population coerced to full panel by caller; focus on the two values
    folds = fold_changes({"CD8": 1.0, "Th1": 1.0, "Th2": 1.0, "Th17": 1.0,
                          "B": 1.0, "NK": 1.0, "Treg": 1.0, "Monocytes": 1.0},
                         {"CD8": 3.0, "Th1": 0.5, "Th2": 1.0, "Th17": 1.0,
                          "B": 1.0, "NK": 1.0, "Treg": 1.0, "Monocytes": 1.0})
    assert folds["CD8"] == 3.0
    assert folds["Th1"] == 0.5


def test_unexpected_population_when_fold_opposes_therapy() -> None:
    folds = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-TNF"]}
    folds["Th17"] = 2.5  # TNF blockade should lower Th17; it rose instead
    assert "Th17" in unexpected_populations(folds, "anti-TNF")
    assert unexpected_populations(folds, "anti-TNF") == ["Th17"]


def test_no_flags_for_normal_response() -> None:
    folds = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}
    folds["CD8"] = 1.5
    assert unexpected_populations(folds, "anti-PD1") == []


def test_evaluate_response_flags_differential_vs_reference() -> None:
    unperturbed = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}

    def with_cd8_fold(fold: float) -> dict:
        perturbed = {pop: 1.0 for pop in TARGET_DIRECTIONS["anti-PD1"]}
        perturbed["CD8"] = fold
        return perturbed

    reference_folds = with_cd8_fold(4.0)
    ref_score = compute_score(reference_folds, "anti-PD1")  # log2(4) = 2.0

    aligned = evaluate_response("anti-PD1", unperturbed, with_cd8_fold(3.8),
                                reference_score=ref_score)
    assert not aligned.is_differential

    divergent = evaluate_response("anti-PD1", unperturbed, with_cd8_fold(0.2),
                                  reference_score=ref_score)
    assert divergent.is_differential
    assert "CD8" in divergent.unexpected