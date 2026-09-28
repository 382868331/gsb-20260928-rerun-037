"""Command line entry point. Sample/reference/parameters come from config."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import __version__
from .config import RunConfig
from .pipeline import run
from .report import (
    build_payload,
    render_text_report,
    write_candidates_tsv,
    write_coverage_tsv,
    write_json,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="variant-pipeline",
        description="Candidate variant summary with evidence records "
        "(research use only).",
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    call = sub.add_parser("call", help="summarize candidate variants")
    call.add_argument("--config", required=True, help="run config JSON path")
    call.add_argument("--reads", help="SAM path (overrides config reads_path)")
    call.add_argument("--outdir", required=True, help="output directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "call":
        cfg = RunConfig.load(args.config)
        result = run(cfg, args.reads)
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        payload = build_payload(cfg, result)
        write_json(payload, outdir / "candidates.json")
        write_candidates_tsv(result, outdir / "candidates.tsv")
        write_coverage_tsv(result, outdir / "coverage.tsv")
        text = render_text_report(payload)
        (outdir / "report.txt").write_text(text, encoding="utf-8")
        print(text)
        print(f"[outputs] {outdir}/candidates.json, candidates.tsv, coverage.tsv")
        return 0
    return 2
