"""Unit tests for the graphe de régulations structure."""

from __future__ import annotations

from bilan_immuno_twin.graph import (
    ACTIVATION_EDGES,
    CYTOKINE_SOURCES,
    CYTOKINES,
    KineticParameters,
    POPULATIONS,
    SUPPRESSION_EDGES,
    activation_edges_for,
    suppression_edges_for,
)


def test_th_subsets_are_separate_nodes() -> None:
    assert {"CD8", "Th1", "Th2", "Th17", "B", "NK", "Treg", "Monocytes"} == set(POPULATIONS)


def test_all_cytokine_nodes_defined_and_quasi_static() -> None:
    assert set(CYTOKINES) == {"IFNg", "IL10", "TNFa", "IL6", "IL4", "IL17"}
    # every cytokine has a producing-population formula
    for cytokine in CYTOKINES:
        assert CYTOKINE_SOURCES[cytokine]
        for source in CYTOKINE_SOURCES[cytokine]:
            assert source in POPULATIONS


def test_edge_sources_and_targets_are_known_nodes() -> None:
    for s, t, k in ACTIVATION_EDGES + SUPPRESSION_EDGES:
        assert s in POPULATIONS + CYTOKINES
        assert t in POPULATIONS
        assert k > 0


def test_classic_cross_regulation_edges_present() -> None:
    # Th1 suppresses Th2, Th2 suppresses Th1
    assert ("Th1", "Th2") in [(s, t) for s, t, _ in SUPPRESSION_EDGES]
    assert ("Th2", "Th1") in [(s, t) for s, t, _ in SUPPRESSION_EDGES]
    # Treg suppresses multiple effector populations
    supp_targets = {t for s, t, _ in SUPPRESSION_EDGES if s == "Treg"}
    assert {"CD8", "Th1", "Th17", "NK", "B"} <= supp_targets


def test_every_population_has_at_least_one_activation_input() -> None:
    act_targets = {t for s, t, _ in ACTIVATION_EDGES}
    assert set(POPULATIONS) <= act_targets


def test_default_parameters_cover_all_edges_and_populations() -> None:
    params = KineticParameters.defaults()
    expected_keys = {(s, t, "act") for s, t, _ in ACTIVATION_EDGES} | {
        (s, t, "sup") for s, t, _ in SUPPRESSION_EDGES
    }
    assert set(params.weights) == expected_keys
    assert set(params.growth) == set(POPULATIONS)
    assert set(params.death) == set(POPULATIONS)
    assert set(params.capacity) == set(POPULATIONS)
    for key, value in params.weights.items():
        assert value > 0
    assert all(v > 0 for v in params.growth.values())
    assert all(v > 0 for v in params.death.values())


def test_activation_and_suppression_helpers() -> None:
    assert activation_edges_for("CD8")  # CD8 has activators (Th1, NK, IFNg)
    for target in POPULATIONS:
        edges = activation_edges_for(target)
        for source, k in edges:
            assert (source, target, k) in ACTIVATION_EDGES