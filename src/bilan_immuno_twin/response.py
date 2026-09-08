"""Score de réponse and differential-response flagging (T4).

The score is deliberately **continuous**, not a binary classifier. For a given
therapy it is a weighted sum of log2 fold-changes of key populations at the end
of a perturbed vs unperturbed simulation:

    score = sum_i target_i * log2(fold_i)

where ``target_i`` is the clinically-desired direction of population ``i`` for
the therapy (+1 = response favors an increase, -1 = response favors a
decrease, 0 = neutral). A positive score is a favorable response, a negative
score an unfavorable one, and near zero means no meaningful response.

A **differential response** flag compares a patient's score against the score
obtained on a reference (typique) bilan with the same therapy: a patient whose
score deviates from the reference beyond a tolerance is flagged as divergent.
Individual populations whose fold direction opposes the therapy's expected
direction, or whose fold magnitude is extreme, are flagged as unexpected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.perturbation import apply_perturbation, get_perturbation
from bilan_immuno_twin.simulation import default_initial_state, simulate

# Desired direction of each population under each therapy.
# +1: response favors an increase, -1: response favors a decrease, 0: neutral.
TARGET_DIRECTIONS: dict[str, dict[str, int]] = {
    "anti-PD1": {"CD8": 1, "Th1": 1, "NK": 1, "B": 0, "Th2": 0,
                 "Th17": 0, "Treg": 0, "Monocytes": 0},
    "anti-TNF": {"Monocytes": -1, "Th17": -1, "Th1": -1, "CD8": 0,
                 "NK": 0, "B": 0, "Th2": 0, "Treg": 0},
    "corticoide": {"Th1": -1, "Th17": -1, "Monocytes": -1, "CD8": 0,
                   "NK": 0, "B": 0, "Th2": 0, "Treg": 0},
}

REFERENCE_BILAN = {
    "CD8": 0.15, "Th1": 0.04, "Th2": 0.06, "Th17": 0.03,
    "B": 0.12, "NK": 0.09, "Treg": 0.02, "Monocytes": 0.11,
}

# A patient is "differential" when their score deviates from the reference score
# by more than this absolute tolerance (in score units).
DIFFERENTIAL_SCORE_TOL = 0.5
# A population is "unexpected" when its fold opposes the desired direction by a
# factor beyond these bounds.
UNEXPECTED_MAX_FOLD = 2.0
UNEXPECTED_MIN_FOLD = 0.5


def _log2(x: float) -> float:
    return 0.0 if x <= 0.0 else math.log2(x)


def fold_changes(unperturbed_values: dict[str, float],
                 perturbed_values: dict[str, float]) -> dict[str, float]:
    """Final fold-change perturbed vs unperturbed per population."""
    folds = {}
    for pop in POPULATIONS:
        base = unperturbed_values[pop]
        treat = perturbed_values[pop]
        if base > 0:
            folds[pop] = treat / base
        else:
            folds[pop] = 0.0
    return folds


def compute_score(folds: dict[str, float], therapy: str) -> float:
    """Continuous score de réponse from per-population fold changes."""
    directions = TARGET_DIRECTIONS[therapy]
    return float(sum(directions[pop] * _log2(fold) for pop, fold in folds.items()))


def score_interpretation(score: float) -> str:
    """Qualitative band for a continuous score (favorable / neutral / unfavorable).

    This complements — but does not replace — the continuous score; it is a
    soft reading aid, not a binary classifier.
    """
    if score > 0.25:
        return "favorable"
    if score < -0.25:
        return "unfavorable"
    return "neutral"


def unexpected_populations(
    folds: dict[str, float],
    therapy: str,
    max_expected_fold: float = UNEXPECTED_MAX_FOLD,
    min_expected_fold: float = UNEXPECTED_MIN_FOLD,
) -> list[str]:
    """Populations whose response diverges from the therapy's expected pattern."""
    flagged = []
    for pop, fold in folds.items():
        desired = TARGET_DIRECTIONS[therapy][pop]
        if desired == 0:
            continue
        if desired == 1 and fold < min_expected_fold:
            flagged.append(pop)
        elif desired == -1 and fold > max_expected_fold:
            flagged.append(pop)
        elif fold > max_expected_fold or fold < min_expected_fold:
            flagged.append(pop)
    return flagged


@dataclass
class Response:
    therapy: str
    score: float
    folds: dict[str, float]
    interpretation: str
    reference_score: float
    is_differential: bool
    unexpected: list[str]


def evaluate_response(
    therapy: str,
    unperturbed_final: dict[str, float],
    perturbed_final: dict[str, float],
    reference_score: float,
    score_tol: float = DIFFERENTIAL_SCORE_TOL,
) -> Response:
    """Bundle a patient's score and flags for a given therapy."""
    folds = fold_changes(unperturbed_final, perturbed_final)
    score = compute_score(folds, therapy)
    return Response(
        therapy=therapy,
        score=score,
        folds=folds,
        interpretation=score_interpretation(score),
        reference_score=reference_score,
        is_differential=abs(score - reference_score) > score_tol,
        unexpected=unexpected_populations(folds, therapy),
    )


def reference_score(therapy: str, horizon: float = 28.0) -> float:
    """Score obtained on the reference bilan with the same therapy.

    Deterministic given (therapy, horizon); used as the "typical expected
    response" the patient's trajectory is compared against.
    """
    params = KineticParameters.defaults()
    x0 = default_initial_state(dict(REFERENCE_BILAN))
    changes = get_perturbation(therapy)
    perturbed = apply_perturbation(params, changes)

    base = simulate(params, x0, horizon=horizon)
    treat = simulate(perturbed, x0, horizon=horizon)
    base_final = {pop: float(base.values[i, -1]) for i, pop in enumerate(POPULATIONS)}
    treat_final = {pop: float(treat.values[i, -1]) for i, pop in enumerate(POPULATIONS)}
    return compute_score(fold_changes(base_final, treat_final), therapy)