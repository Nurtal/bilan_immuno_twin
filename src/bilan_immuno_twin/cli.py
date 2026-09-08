"""CLI entry point for the bilan immuno twin (ADR-0004: CLI only)."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Sequence

import numpy as np

from bilan_immuno_twin.bilan import BilanError, parse_bilan
from bilan_immuno_twin.graph import KineticParameters, POPULATIONS
from bilan_immuno_twin.perturbation import PerturbationError, apply_perturbation, get_perturbation
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


def _trajectory_to_dict(times, values) -> dict:
    populations = {}
    for i, pop in enumerate(POPULATIONS):
        populations[pop] = {"label": POPULATION_LABELS[pop], "times": times.tolist(), "values": values[i].tolist()}
    return {
        "populations": populations,
        "t_units": "days",
        "stopped_early": False,
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
    p_sim.set_defaults(func=_cmd_simulate)
    return parser


def _cmd_bilan(args: argparse.Namespace) -> int:
    populations = parse_bilan(args.bilan_file)
    print(json.dumps({"populations": populations}, sort_keys=True))
    return 0


def _cmd_simulate(args: argparse.Namespace) -> int:
    bilan = parse_bilan(args.bilan_file)
    x0 = default_initial_state(bilan)
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
            baseline, x0, horizon=args.horizon, adaptive=args.adaptive
        )
        output["perturbed"] = _simulate_output(
            perturbed, x0, horizon=args.horizon, adaptive=args.adaptive
        )
        output["comparison"] = _comparison(output["unperturbed"], output["perturbed"])
    else:
        output = _simulate_output(baseline, x0, horizon=args.horizon, adaptive=args.adaptive)

    print(json.dumps(output))
    return 0


def _simulate_output(params: KineticParameters, x0: np.ndarray, horizon: float,
                     adaptive: bool) -> dict:
    trj = simulate(params, x0, horizon=horizon, adaptive=adaptive)
    output = _trajectory_to_dict(trj.times, trj.values)
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