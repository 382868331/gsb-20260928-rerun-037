"""Command line entry point. Sample/reference/parameters come from config."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from . import __version__
from .config import RunConfig
from .diff import compare_versions, render_diff_text
from .pipeline import run
from .report import (
    build_payload,
    render_text_report,
    write_candidates_tsv,
    write_coverage_tsv,
    write_json,
)
from .store import OUTPUT_FILENAMES, RunStore, StoreError, VerificationError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="variant-pipeline",
        description="Candidate variant summary with evidence records "
        "(research use only).",
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    call = sub.add_parser(
        "call", help="analyse and publish a new immutable result version")
    call.add_argument("--config", required=True, help="run config JSON path")
    call.add_argument("--reads", help="SAM path (overrides config reads_path)")
    call.add_argument("--outdir", required=True,
                      help="versioned run store root directory")
    call.add_argument("--force", action="store_true",
                      help="publish a new version even if an identical "
                      "request was already published")

    list_cmd = sub.add_parser("list", help="list published versions")
    list_cmd.add_argument("--outdir", required=True)

    compare = sub.add_parser(
        "compare", help="diff candidates between two versions")
    compare.add_argument("--outdir", required=True)
    compare.add_argument("--from", dest="from_version", required=True)
    compare.add_argument("--to", dest="to_version", required=True)
    compare.add_argument("--json-out",
                         help="also write the structured diff JSON here")

    export = sub.add_parser(
        "export", help="export a historical version after integrity check")
    export.add_argument("--outdir", required=True)
    export.add_argument("--version", required=True)
    export.add_argument("--dest", required=True,
                        help="empty/non-existent destination directory")

    verify = sub.add_parser(
        "verify", help="verify stored versions were not rewritten")
    verify.add_argument("--outdir", required=True)
    verify.add_argument("--version",
                        help="verify a single version (default: all)")
    return parser


def _cmd_call(args) -> int:
    cfg = RunConfig.load(args.config)
    reads = args.reads or cfg.reads_path
    if not reads:
        raise StoreError(
            "no reads path: set reads_path in config or pass --reads")
    store = RunStore(args.outdir)
    config_dict = cfg.to_dict()
    if args.reads:
        config_dict["reads_path"] = args.reads
    fingerprint = store.fingerprint(
        config_dict=config_dict,
        reads_path=reads,
        fasta_path=cfg.reference.fasta_path,
    )
    if not args.force:
        duplicate = store.find_duplicate(fingerprint)
        if duplicate is not None:
            print(f"[duplicate] 相同输入与规则的请求已发布为 "
                  f"{duplicate.version}，不重复分析；"
                  f"结果见 {duplicate.path}")
            print(f"          如需强制生成新版本，请加 --force")
            return 0

    result = run(cfg, args.reads)
    payload = build_payload(cfg, result)
    text = render_text_report(payload)
    with tempfile.TemporaryDirectory(
            prefix=".call-", dir=str(Path(args.outdir).parent)
            if Path(args.outdir).parent.exists() else None) as tmp:
        tmpdir = Path(tmp)
        write_json(payload, tmpdir / "candidates.json")
        write_candidates_tsv(result, tmpdir / "candidates.tsv")
        write_coverage_tsv(result, tmpdir / "coverage.tsv")
        (tmpdir / "report.txt").write_text(text, encoding="utf-8")
        outputs = {name: tmpdir / name for name in OUTPUT_FILENAMES}
        info = store.publish(
            outputs,
            config_dict=config_dict,
            reads_path=reads,
            fasta_path=cfg.reference.fasta_path,
            fingerprint=fingerprint,
        )
    print(text)
    print(f"[version] 已发布 {info.version} -> {info.path}")
    print(f"[outputs] {info.path}/candidates.json, candidates.tsv, "
          f"coverage.tsv, report.txt")
    return 0


def _cmd_list(args) -> int:
    store = RunStore(args.outdir)
    infos = store.list_versions()
    if not infos:
        print("无已发布版本")
        return 0
    for info in infos:
        manifest = info.manifest
        sample_id = (
            manifest.get("inputs", {}).get("config", {})
            .get("sample", {}).get("sample_id", "-")
        )
        print(f"{info.version}	{manifest.get('created_at', '-')}"
              f"	{manifest.get('analysis_version', '-')}"
              f"	sample={sample_id}")
    return 0
def _cmd_compare(args) -> int:
    store = RunStore(args.outdir)
    diff = compare_versions(store, args.from_version, args.to_version)
    print(render_diff_text(diff))
    if args.json_out:
        Path(args.json_out).write_text(
            json.dumps(diff, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")
        print(f"[diff] {args.json_out}")
    return 0


def _cmd_export(args) -> int:
    store = RunStore(args.outdir)
    dest = store.export(args.version, args.dest)
    print(f"[export] {args.version} 完整性核对通过，已导出到 {dest}")
    return 0


def _cmd_verify(args) -> int:
    store = RunStore(args.outdir)
    if args.version:
        store.verify(args.version)
        print(f"[verify] {args.version} 内容完整，未被改写")
        return 0
    infos = store.verify_all()
    if not infos:
        print("无已发布版本")
        return 0
    for info in infos:
        print(f"[verify] {info.version} 内容完整，未被改写")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handlers = {
        "call": _cmd_call,
        "list": _cmd_list,
        "compare": _cmd_compare,
        "export": _cmd_export,
        "verify": _cmd_verify,
    }
    try:
        return handlers[args.command](args)
    except (StoreError, VerificationError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

