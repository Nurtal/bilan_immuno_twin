"""Graphe de régulations: fixed topology, quasi-static cytokines, default parameters.

The graph has two kinds of nodes:
  - **Populations** (variables d'état): CD8, Th1, Th2, Th17, B, NK, Treg, Mono.
    Each is a separate node of the ODE system.
  - **Cytokines** (variables quasi-statiques, ADR-0005): computed algebraically from
    the producing populations at each timestep, with no ODE of their own.

Directed edges encode regulatory interactions from immunology literature:
activation and suppression. Topology is fixed; edge weights are patient-specific
kinetic parameters.
"""

from __future__ import annotations

from dataclasses import dataclass, field

POPULATIONS = ["CD8", "Th1", "Th2", "Th17", "B", "NK", "Treg", "Monocytes"]

CYTOKINES = ["IFNg", "IL10", "TNFa", "IL6", "IL4", "IL17"]

ALL_NODES = POPULATIONS + CYTOKINES

HILL_N = 2.0

# --- Fixed topology ----------------------------------------------------------

# Activation edges: (source, target, K) — Hill threshold in fraction units.
# Literature grounding: Th1/IFN-γ help CD8+ cytotoxics; NK co-stimulates CD8+;
# IL-12 from monocytes/macrophages drives Th1; IL-4 from B cells drives Th2;
# IL-6/TGF-β from monocytes drives Th17 and potentiates inflammation (TNF-α);
# helper subsets drive B cells and NK; IFN-γ and IL-4 feed back on their lineage.
ACTIVATION_EDGES = [
    ("Th1", "CD8", 0.12),
    ("NK", "CD8", 0.15),
    ("IFNg", "CD8", 0.18),
    ("CD8", "Th1", 0.10),
    ("Monocytes", "Th1", 0.15),
    ("IFNg", "Th1", 0.20),
    ("B", "Th2", 0.12),
    ("IL4", "Th2", 0.15),
    ("Monocytes", "Th17", 0.15),
    ("IL6", "Th17", 0.18),
    ("Th2", "B", 0.12),
    ("IL4", "B", 0.15),
    ("Th1", "NK", 0.12),
    ("Monocytes", "NK", 0.12),
    ("IFNg", "NK", 0.18),
    ("Th2", "Treg", 0.12),
    ("IL10", "Treg", 0.15),
    ("Th1", "Monocytes", 0.12),
    ("Th17", "Monocytes", 0.12),
    ("TNFa", "Monocytes", 0.15),
    ("IL6", "Monocytes", 0.15),
]

# Suppression edges: (source, target, K).
# Literature grounding: Tregs suppress effector responses (CD8, Th1, Th17, NK, B)
# via IL-10/contact; Th1 (IFN-γ) cross-suppresses Th2; Th2 (IL-4/IL-10)
# cross-suppresses Th1; Th17 curbs Treg induction.
SUPPRESSION_EDGES = [
    ("Treg", "CD8", 0.08),
    ("Treg", "Th1", 0.08),
    ("Treg", "Th17", 0.08),
    ("Treg", "NK", 0.08),
    ("Treg", "B", 0.08),
    ("Th1", "Th2", 0.15),
    ("IL10", "Th1", 0.15),
    ("Th2", "Th1", 0.20),
    ("Th17", "Treg", 0.18),
]

# Quasi-static cytokine production (ADR-0005): each cytokine is a linear function
# of its producing populations. Coefficients in shared arbitrary units.
CYTOKINE_SOURCES = {
    "IFNg": {"Th1": 1.0, "NK": 0.5, "CD8": 0.3},
    "IL10": {"Th2": 1.0, "Treg": 1.2},
    "TNFa": {"Th1": 1.0, "Th17": 0.8, "Monocytes": 1.0},
    "IL6": {"Monocytes": 1.0, "Th2": 0.3},
    "IL4": {"Th2": 1.0},
    "IL17": {"Th17": 1.0},
}

# Edge identifiers: kind in {"act", "sup"}.
EDGE_KIND: dict[tuple[str, str], str] = {}
for _s, _t, _k in ACTIVATION_EDGES:
    EDGE_KIND[(_s, _t)] = "act"
for _s, _t, _k in SUPPRESSION_EDGES:
    EDGE_KIND[(_s, _t)] = "sup"


@dataclass
class KineticParameters:
    """Patient-specific kinetic parameters of the graphe de régulations.

    - growth: basal proliferation rate per population (per day)
    - death:  basal death rate per population (per day)
    - capacity: carrying capacity per population (fraction of total immune cells)
    - basal: small basal activation feeding each population's growth
    - weights: edge weights keyed by (source, target, kind); the patient-specific
      reactivity of each regulatory interaction
    """

    growth: dict[str, float] = field(default_factory=dict)
    death: dict[str, float] = field(default_factory=dict)
    capacity: dict[str, float] = field(default_factory=dict)
    basal: dict[str, float] = field(default_factory=dict)
    weights: dict[tuple[str, str, str], float] = field(default_factory=dict)

    def edge_weight(self, source: str, target: str) -> float:
        return self.weights[(source, target, EDGE_KIND[(source, target)])]

    @classmethod
    def defaults(cls) -> "KineticParameters":
        growth = {p: 0.35 for p in POPULATIONS}
        death = {p: 0.08 for p in POPULATIONS}
        capacity = {p: 1.0 for p in POPULATIONS}
        basal = {p: 0.12 for p in POPULATIONS}
        weights = {}
        for s, t, _k in ACTIVATION_EDGES:
            weights[(s, t, "act")] = 0.6
        for s, t, _k in SUPPRESSION_EDGES:
            weights[(s, t, "sup")] = 0.6
        return cls(growth=growth, death=death, capacity=capacity,
                   basal=basal, weights=weights)


def activation_edges_for(target: str) -> list[tuple[str, float]]:
    return [(s, k) for s, t, k in ACTIVATION_EDGES if t == target]


def suppression_edges_for(target: str) -> list[tuple[str, float]]:
    return [(s, k) for s, t, k in SUPPRESSION_EDGES if t == target]