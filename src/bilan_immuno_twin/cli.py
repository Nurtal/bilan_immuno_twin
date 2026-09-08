"""CLI entry point for the bilan immuno twin (ADR-0004: CLI only)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

from bilan_immuno_twin.bilan import BilanError, parse_bilan
from bilan_immuno_twin.calibration import calibrate, parameters_from_growth
from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.perturbation import PerturbationError, apply_perturbation, get_perturbation
from bilan_immuno_twin.response import evaluate_response, reference_score
from bilan_immuno_twin.simulation import default_initial_state, simulate

POPULATION_LABELS = {
    "CD8": "CD8+ T",
    "Th1": "Th1",
    "Th2": "Th2",
    "Th17": "Th17",
    "B": "B",
    "NK": "NK",
    "Treg": "Treg",
    "Monocytes": "Monocytes",
}


def _trajectory_to_dict(times, values, ci: dict | None = None) -> dict:
    populations = {}
    for i, pop in enumerate(POPULATIONS):
        entry = {"label": POPULATION_LABELS[pop], "times": times.tolist(), "values": values[i].tolist()}
        if ci is not None:
            entry["ci_lo"] = ci[pop]["ci_lo"]
            entry["ci_hi"] = ci[pop]["ci_hi"]
        populations[pop] = entry
    return {
        "populations": populations,
        "t_units": "days",
        "stopped_early": False,
    }


def _bootstrap_ci_band(
    bootstrap_values: list,
    x0: np.ndarray,
    horizon: float,
    perturb: str | None,
    max_draws: int = 30,
) -> dict:
    """95% percentile band over the bootstrap parameter draws for each population.

    Each bootstrap growth draw is simulated over the same horizon (with the same
    perturbation applied) and interpolated onto a shared time grid; the 2.5th
    and 97.5th percentiles across draws form the band.
    """
    grid = np.linspace(0.0, horizon, 200)
    draws = bootstrap_values[:max_draws]
    per_pop = {pop: np.empty((len(draws), grid.size), dtype=float) for pop in POPULATIONS}
    for d, draw in enumerate(draws):
        params = parameters_from_growth(dict(zip(POPULATIONS, draw)))
        if perturb:
            params = apply_perturbation(params, get_perturbation(perturb))
        trj = simulate(params, x0, horizon=horizon)
        for i, pop in enumerate(POPULATIONS):
            per_pop[pop][d] = np.interp(grid, trj.times, trj.values[i])
    return {
        pop: {
            "ci_lo": np.percentile(per_pop[pop], 2.5, axis=0).tolist(),
            "ci_hi": np.percentile(per_pop[pop], 97.5, axis=0).tolist(),
        }
        for pop in POPULATIONS
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bilan",
        description="Jumeau numérique immunitaire — digital twin of the immune system",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_bilan = sub.add_parser("bilan", help="read a bilan immunologique and emit the measured populations")
    p_bilan.add_argument("bilan_file", help="path to the bilan immunologique file (JSON)")
    p_bilan.set_defaults(func=_cmd_bilan)

    p_sim = sub.add_parser("simulate", help="simulate immune population dynamics from a bilan")
    p_sim.add_argument("bilan_file", help="path to the bilan immunologique file (JSON)")
    p_sim.add_argument("--horizon", type=float, default=28.0, help="simulation horizon in days")
    p_sim.add_argument("--adaptive", action="store_true", help="stop early once a steady state is reached")
    p_sim.add_argument("--perturb", metavar="THERAPY",
                       help="apply an immunotherapy perturbation (anti-PD1, anti-TNF, corticoide)")
    p_sim.add_argument("--params", metavar="JSON",
                       help="calibrated parameters file (calibrate output or {\"growth\": {...}})")
    p_sim.add_argument("--with-ci", action="store_true",
                       help="emit 95% confidence bands from the bootstrap draws in --params")
    p_sim.set_defaults(func=_cmd_simulate)

    p_cal = sub.add_parser("calibrate", help="MAP-calibrate patient-specific kinetic parameters from a bilan")
    p_cal.add_argument("bilan_file", help="path to the bilan immunologique file (JSON)")
    p_cal.add_argument("--bootstrap", type=int, default=100, help="number of bootstrap resamples")
    p_cal.add_argument("--seed", type=int, default=None, help="random seed for reproducibility")
    p_cal.add_argument("--out", metavar="JSON", default=None,
                       help="write the calibration result to a JSON file")
    p_cal.set_defaults(func=_cmd_calibrate)
    return parser


def _cmd_bilan(args: argparse.Namespace) -> int:
    populations = parse_bilan(args.bilan_file)
    print(json.dumps({"populations": populations}, sort_keys=True))
    return 0


def _calibration_params_from_file(path: str) -> tuple[KineticParameters, dict | None]:
    """Load calibrated parameters (and optional bootstrap draws) from a JSON file.

    The file is either a `calibrate` CLI output or a plain {"growth": {...}} map.
    """
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    growth = raw.get("growth")
    if growth is None and isinstance(raw.get("map"), dict):
        growth = raw["map"].get("growth")
    if not isinstance(growth, dict):
        raise BilanError(f"params file {path}: expected a 'growth' map of calibrated parameters")
    bootstrap_values = raw.get("bootstrap_values")
    return parameters_from_growth(growth), bootstrap_values


def _cmd_simulate(args: argparse.Namespace) -> int:
    bilan = parse_bilan(args.bilan_file)
    x0 = default_initial_state(bilan)

    bootstrap_values = None
    if args.params:
        baseline, bootstrap_values = _calibration_params_from_file(args.params)
    else:
        baseline = KineticParameters.defaults()

    output = {}
    if args.perturb:
        changes = get_perturbation(args.perturb)
        perturbed = apply_perturbation(baseline, changes)
        output["perturbation"] = args.perturb
        output["parameter_changes"] = {
            identifier: {"delta": delta, "interpretation": f"{delta * 100:+.0f}%"}
            for identifier, delta in changes.items()
        }
        output["unperturbed"] = _simulate_output(
            baseline, x0, horizon=args.horizon, adaptive=args.adaptive,
            bootstrap_values=bootstrap_values if args.with_ci else None,
            perturb=args.perturb,
        )
        output["perturbed"] = _simulate_output(
            perturbed, x0, horizon=args.horizon, adaptive=args.adaptive,
            bootstrap_values=bootstrap_values if args.with_ci else None,
            perturb=args.perturb,
        )
        output["comparison"] = _comparison(output["unperturbed"], output["perturbed"])
        output["response"] = _response_summary(args.perturb, output, horizon=args.horizon)
    else:
        output = _simulate_output(
            baseline, x0, horizon=args.horizon, adaptive=args.adaptive,
            bootstrap_values=bootstrap_values if args.with_ci else None,
            perturb=args.perturb,
        )

    print(json.dumps(output))
    return 0


def _cmd_calibrate(args: argparse.Namespace) -> int:
    bilan = parse_bilan(args.bilan_file)
    result = calibrate(bilan, n_bootstrap=args.bootstrap, seed=args.seed)
    output = {
        "map": {
            "growth": {name: float(value) for name, value in result.growth_estimates.items()},
        },
        "confidence_intervals": {
            name: [float(lo), float(hi)]
            for name, (lo, hi) in result.confidence_intervals.items()
        },
        "bootstrap_values": result.bootstrap_values.tolist(),
        "n_bootstrap": args.bootstrap,
        "converged": result.converged,
        "mse": result.mse,
    }
    if args.out:
        Path(args.out).write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output))
    return 0


def _simulate_output(params: KineticParameters, x0: np.ndarray, horizon: float,
                     adaptive: bool, bootstrap_values: list | None = None,
                     perturb: str | None = None) -> dict:
    trj = simulate(params, x0, horizon=horizon, adaptive=adaptive)
    ci = None
    if bootstrap_values:
        ci = _bootstrap_ci_band(bootstrap_values, x0, horizon, perturb)
    output = _trajectory_to_dict(trj.times, trj.values, ci=ci)
    output["stopped_early"] = trj.stopped_early
    return output


def _comparison(unperturbed: dict, perturbed: dict) -> dict:
    """Fold-change of final population values perturbed vs unperturbed."""
    comparison = {}
    for pop in POPULATIONS:
        base = unperturbed["populations"][pop]["values"][-1]
        treat = perturbed["populations"][pop]["values"][-1]
        fold = treat / base if base > 0 else None
        direction = "up" if (fold or 1.0) > 1.0 else ("down" if (fold or 1.0) < 1.0 else "flat")
        comparison[pop] = {"fold_change": fold, "direction": direction}
    return comparison


def _response_summary(therapy: str, output: dict, horizon: float) -> dict:
    unperturbed = output["unperturbed"]["populations"]
    perturbed = output["perturbed"]["populations"]
    base_final = {pop: unperturbed[pop]["values"][-1] for pop in POPULATIONS}
    treat_final = {pop: perturbed[pop]["values"][-1] for pop in POPULATIONS}
    response = evaluate_response(
        therapy,
        base_final,
        treat_final,
        reference_score=reference_score(therapy, horizon=horizon),
    )
    return {
        "score": response.score,
        "therapy": therapy,
        "interpretation": response.interpretation,
        "reference_score": response.reference_score,
        "is_differential": response.is_differential,
        "unexpected_populations": response.unexpected,
        "folds": response.folds,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (BilanError, PerturbationError, ValueError) as exc:
        print(f"bilan: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())