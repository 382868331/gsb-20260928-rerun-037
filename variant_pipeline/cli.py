"""Command line entry point.

The call subcommand publishes every run as a new, immutable versioned
batch under <outdir>/batches/ and never overwrites previous batches. The
list/diff/export/verify subcommands operate only on frozen historical batch
files; they never re-run analysis on old inputs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import __version__
from . import store
from .config import RunConfig
from .diff import diff_batches, render_diff_text
from .pipeline import run
from .report import (
    build_payload,
    candidates_tsv_text,
    coverage_tsv_text,
    render_text_report,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="variant-pipeline",
        description="Versioned candidate variant summary with evidence "
        "records (research use only).",
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    call = sub.add_parser("call", help="analyze inputs and publish a new batch")
    call.add_argument("--config", required=True, help="run config JSON path")
    call.add_argument("--reads", help="SAM path (overrides config reads_path)")
    call.add_argument(
        "--outdir",
        required=True,
        help="versioned output root; the batch is stored under batches/<id>",
    )

    lst = sub.add_parser("list", help="list historical batches")
    lst.add_argument("--outdir", required=True, help="versioned output root")
    lst.add_argument("--json", action="store_true", help="emit JSON")

    dif = sub.add_parser("diff", help="compare two historical batches")
    dif.add_argument("--outdir", required=True, help="versioned output root")
    dif.add_argument("--from", dest="from_batch", default="latest",
                     help="old batch id/prefix (default: latest)")
    dif.add_argument("--to", dest="to_batch", required=True,
                     help="new batch id/prefix")
    dif.add_argument("--json", action="store_true", help="emit JSON")

    exp = sub.add_parser("export", help="export a historical batch as-is")
    exp.add_argument("--outdir", required=True, help="versioned output root")
    exp.add_argument("--batch", default="latest", help="batch id/prefix")
    exp.add_argument("--dest", required=True, help="export destination directory")
    exp.add_argument("--overwrite", action="store_true")

    ver = sub.add_parser(
        "verify", help="verify stored hashes without re-analyzing inputs"
    )
    ver.add_argument("--outdir", required=True, help="versioned output root")
    ver.add_argument("--batch", default="latest", help="batch id/prefix")
    return parser


def _cmd_call(args) -> int:
    cfg = RunConfig.load(args.config)
    inputs = store.load_input_set(cfg, args.config, args.reads)
    result = run(cfg, args.reads)

    root = Path(args.outdir)
    index = store.load_index(root)
    duplicates = store.find_duplicate_batches(root, inputs.content_id())
    payload = build_payload(
        cfg,
        result,
        batch={
            "content_id": inputs.content_id(),
            "input_hashes": inputs.input_hashes,
            "input_sources": inputs.input_sources,
            "duplicate_of": duplicates,
        },
    )
    payload_text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    report_text = render_text_report(payload)
    artifacts = {
        "candidates.json": payload_text.encode("utf-8"),
        "candidates.tsv": candidates_tsv_text(result).encode("utf-8"),
        "coverage.tsv": coverage_tsv_text(result).encode("utf-8"),
        "report.txt": report_text.encode("utf-8"),
    }
    published = store.publish_batch(root, cfg, inputs, payload, artifacts)
    batch_id = published["manifest"]["batch_id"]
    print(report_text, end="")
    if published["manifest"]["duplicate_of"]:
        print(
            "[duplicate] identical analysis_version+inputs already published: "
            + ", ".join(published["manifest"]["duplicate_of"])
            + "; previous batches left untouched"
        )
    print(f"[batch] {batch_id}")
    print(
        f"[outputs] {published['batch_dir']}/"
        "candidates.json, candidates.tsv, coverage.tsv, report.txt"
    )
    return 0


def _cmd_list(args) -> int:
    index = store.load_index(args.outdir)
    batches = index.get("batches", [])
    if args.json:
        print(json.dumps(index, indent=2, ensure_ascii=False))
        return 0
    if not batches:
        print("(no batches)")
        return 0
    for b in batches:
        dup = " DUPLICATE" if b.get("duplicate_of") else ""
        print(
            f"{b['batch_id']}  {b['analysis_version']}  "
            f"{b['sample_id']}  candidates={b['n_candidates']}"
            f"  assembly={b['reference_assembly']}{dup}"
        )
    return 0


def _cmd_diff(args) -> int:
    result = diff_batches(args.outdir, args.from_batch, args.to_batch)
    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        print(render_diff_text(result), end="")
    return 0


def _cmd_export(args) -> int:
    batch_id = store.resolve_batch_id(args.outdir, args.batch)
    dest = store.export_batch(args.outdir, batch_id, args.dest, args.overwrite)
    print(f"[exported] {batch_id} -> {dest} (stored copy, not re-analyzed)")
    return 0


def _cmd_verify(args) -> int:
    batch_id = store.resolve_batch_id(args.outdir, args.batch)
    result = store.verify_batch(args.outdir, batch_id)
    for rel, check in result["files"].items():
        status = (
            "OK"
            if check["ok"]
            else f"FAIL ({check.get('reason', 'hash mismatch')})"
        )
        print(f"  {status}  {rel}")
    print(f"[verify] {batch_id}: {'OK' if result['ok'] else 'FAILED'}")
    return 0 if result["ok"] else 1


_HANDLERS = {
    "call": _cmd_call,
    "list": _cmd_list,
    "diff": _cmd_diff,
    "export": _cmd_export,
    "verify": _cmd_verify,
}


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return _HANDLERS[args.command](args)
