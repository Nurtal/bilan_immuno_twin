"""CLI entry point for the bilan immuno twin (ADR-0004: CLI only)."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Sequence

from bilan_immuno_twin.bilan import BilanError, parse_bilan


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bilan",
        description="Jumeau numérique immunitaire — digital twin of the immune system",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_bilan = sub.add_parser("bilan", help="read a bilan immunologique and emit the measured populations")
    p_bilan.add_argument("bilan_file", help="path to the bilan immunologique file (JSON)")
    p_bilan.set_defaults(func=_cmd_bilan)
    return parser


def _cmd_bilan(args: argparse.Namespace) -> int:
    populations = parse_bilan(args.bilan_file)
    print(json.dumps({"populations": populations}, sort_keys=True))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except BilanError as exc:
        print(f"bilan: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())