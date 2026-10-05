"""Command line entry point. Sample/reference/parameters come from config."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from . import __version__
from .combination import (
    adopt_combination,
    load_combination_constraints,
    search_combinations,
)
from .config import RunConfig
from .diff import compare_versions, render_diff_text
from .interpretation import (
    load_selection,
    save_selection,
    selected_view,
    withdraw_selection,
)
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
    call.add_argument("--extra-reads", action="append", default=[],
                      help="additional SAM batch; repeatable")
    call.add_argument("--associations", help="JSON file of explicit same-molecule fragment relations")
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

    select = sub.add_parser(
        "select", help="save a provisional candidate interpretation")
    select.add_argument("--outdir", required=True)
    select.add_argument("--version", required=True)
    select.add_argument("--variant", required=True)
    select.add_argument("--evidence-unit", help="association unit ID adopted for this interpretation")
    select.add_argument("--relation-id")
    select.add_argument("--rationale")

    withdraw = sub.add_parser(
        "withdraw", help="withdraw a provisional interpretation choice")
    withdraw.add_argument("--outdir", required=True)
    withdraw.add_argument("--version", required=True)

    view = sub.add_parser(
        "selected-view", help="show evidence under a provisional choice")
    view.add_argument("--outdir", required=True)
    view.add_argument("--version", required=True)
    view.add_argument("--json-out")

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
    combos = sub.add_parser(
        "combinations",
        help="search compatible relation combinations under declared "
        "co-required / mutually-exclusive constraints (does not confirm "
        "any relation)")
    combos.add_argument("--outdir", required=True)
    combos.add_argument("--version", required=True)
    combos.add_argument("--constraints", required=True,
                        help="JSON file of combination constraints")
    combos.add_argument("--json-out")

    adopt = sub.add_parser(
        "adopt-combination",
        help="explicitly adopt one candidate combination plan; records "
        "the selected sources without modifying the base version")
    adopt.add_argument("--outdir", required=True)
    adopt.add_argument("--version", required=True)
    adopt.add_argument("--constraints", required=True)
    adopt.add_argument("--plan", required=True, help="candidate plan ID")
    adopt.add_argument("--rationale")
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
    if args.associations:
        config_dict["associations_path"] = args.associations
    fingerprint = store.fingerprint(
        config_dict=config_dict,
        reads_path=reads,
        fasta_path=cfg.reference.fasta_path,
        extra_read_paths=args.extra_reads,
        associations_path=args.associations or cfg.associations_path,
    )
    if not args.force:
        duplicate = store.find_duplicate(fingerprint)
        if duplicate is not None:
            print(f"[duplicate] 相同输入与规则的请求已发布为 "
                  f"{duplicate.version}，不重复分析；"
                  f"结果见 {duplicate.path}")
            print(f"          如需强制生成新版本，请加 --force")
            return 0

    read_paths = [reads] + args.extra_reads
    if len(read_paths) == 1:
        result = run(cfg, reads, args.associations or cfg.associations_path)
    else:
        from .pipeline import run_many
        result = run_many(cfg, read_paths, args.associations or cfg.associations_path)
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
            extra_read_paths=args.extra_reads,
            associations_path=args.associations or cfg.associations_path,
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


def _cmd_select(args) -> int:
    selection = save_selection(
        RunStore(args.outdir), args.version, args.variant,
        evidence_unit_id=args.evidence_unit, relation_id=args.relation_id,
        rationale=args.rationale)
    print(f"[selection] 已临时选择 {selection.selected_variant_id}; "
          f"原分析版本 {selection.base_version} 未修改")
    return 0


def _cmd_withdraw(args) -> int:
    selection = withdraw_selection(RunStore(args.outdir), args.version)
    print(f"[withdraw] 已撤回 {selection.base_version} 的临时选择，原版本未修改")
    return 0


def _cmd_selected_view(args) -> int:
    store = RunStore(args.outdir)
    decision = load_selection(store, args.version)
    if decision is None:
        raise StoreError(f"no provisional selection for {args.version}")
    payload = selected_view(
        store.verify(args.version),
        decision["selected_variant_id"], decision)
    text = json.dumps(payload, indent=2, ensure_ascii=False)
    print(text)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n", encoding="utf-8")
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


def _cmd_combinations(args) -> int:
    store = RunStore(args.outdir)
    info = store.verify(args.version)
    payload = json.loads(
        (info.path / "candidates.json").read_text(encoding="utf-8"))
    constraints = load_combination_constraints(args.constraints)
    result = search_combinations(
        payload.get("evidence_units", []), constraints)
    text = json.dumps(result.to_dict(), indent=2, ensure_ascii=False)
    print(text)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n", encoding="utf-8")
    if result.status != "ok":
        print("[combinations] 约束互相冲突，无可采用组合；"
              "冲突约束与端点见 conflicts", file=sys.stderr)
        return 1
    print("[combinations] 候选方案仅供比较，未确认任何关系；"
          "明确采用请用 adopt-combination --plan PLAN_ID")
    return 0


def _cmd_adopt_combination(args) -> int:
    store = RunStore(args.outdir)
    constraints = load_combination_constraints(args.constraints)
    record = adopt_combination(
        store, args.version, constraints, args.plan,
        rationale=args.rationale)
    print(f"[combination] 已采用 {record['plan_id']} "
          f"(关系 {len(record['association_ids'])} 条, "
          f"来源 {record['independent_source_count']} 个: "
          f"{', '.join(record['selected_sources'])}); "
          f"原分析版本 {record['base_version']} 未修改")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handlers = {
        "call": _cmd_call,
        "list": _cmd_list,
        "compare": _cmd_compare,
        "select": _cmd_select,
        "withdraw": _cmd_withdraw,
        "selected-view": _cmd_selected_view,
        "export": _cmd_export,
        "verify": _cmd_verify,
        "combinations": _cmd_combinations,
        "adopt-combination": _cmd_adopt_combination,
    }
    try:
        return handlers[args.command](args)
    except (StoreError, VerificationError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
