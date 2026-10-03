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
from .interpretation import (
    link_view,
    list_link_selections,
    load_link_selection,
    load_selection,
    save_link_selection,
    save_selection,
    selected_view,
    withdraw_link_selection,
    withdraw_selection,
)
from .linkage import load_links, recount_payload
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
    call.add_argument("--links", action="append", default=[],
                      help="linkage JSON declaring fragments from the same "
                      "original molecule; repeatable")
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

    link_view_cmd = sub.add_parser(
        "link-view", help="show one linkage unit and its adoption state")
    link_view_cmd.add_argument("--outdir", required=True)
    link_view_cmd.add_argument("--version", required=True)
    link_view_cmd.add_argument("--link", required=True)
    link_view_cmd.add_argument("--json-out")

    link_select = sub.add_parser(
        "link-select",
        help="provisionally adopt one interpretation for a linkage unit")
    link_select.add_argument("--outdir", required=True)
    link_select.add_argument("--version", required=True)
    link_select.add_argument("--link", required=True)
    link_select.add_argument("--variant", required=True)
    link_select.add_argument("--rationale", required=True)

    link_withdraw = sub.add_parser(
        "link-withdraw", help="withdraw a linkage-unit adoption")
    link_withdraw.add_argument("--outdir", required=True)
    link_withdraw.add_argument("--version", required=True)
    link_withdraw.add_argument("--link", required=True)

    link_compare = sub.add_parser(
        "link-compare",
        help="compare molecule-level support under an alternative "
        "linkage rule without modifying the base version")
    link_compare.add_argument("--outdir", required=True)
    link_compare.add_argument("--version", required=True)
    link_compare.add_argument("--links", required=True,
                              help="alternative linkage JSON file")
    link_compare.add_argument("--json-out")
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
    link_specs = []
    for link_path in args.links:
        link_specs.extend(load_links(link_path))
    fingerprint = store.fingerprint(
        config_dict=config_dict,
        reads_path=reads,
        fasta_path=cfg.reference.fasta_path,
        extra_read_paths=args.extra_reads,
        link_paths=args.links,
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
        result = run(cfg, reads, links=link_specs)
    else:
        from .pipeline import run_many
        result = run_many(cfg, read_paths, links=link_specs)
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
            link_paths=args.links,
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
    selection = save_selection(RunStore(args.outdir), args.version, args.variant)
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
    payload = selected_view(store.verify(args.version),
                            decision["selected_variant_id"])
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


def _cmd_link_view(args) -> int:
    store = RunStore(args.outdir)
    info = store.verify(args.version)
    adoption = load_link_selection(store, args.version, args.link)
    payload = link_view(info, args.link, adoption=adoption)
    text = json.dumps(payload, indent=2, ensure_ascii=False)
    print(text)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n", encoding="utf-8")
    return 0


def _cmd_link_select(args) -> int:
    store = RunStore(args.outdir)
    record = save_link_selection(
        store, args.version, args.link, args.variant, args.rationale)
    print(f"[link-select] 已临时采用 {record['selected_variant_id']} "
          f"解释关联单元 {record['link_id']}; "
          f"排除候选={record['excluded_variant_ids']}; "
          f"原分析版本 {record['base_version']} 未修改")
    return 0


def _cmd_link_withdraw(args) -> int:
    store = RunStore(args.outdir)
    record = withdraw_link_selection(store, args.version, args.link)
    print(f"[link-withdraw] 已撤回 {record['base_version']} 关联单元 "
          f"{record['link_id']} 的临时采用，恢复原有歧义视图")
    return 0


def _cmd_link_compare(args) -> int:
    store = RunStore(args.outdir)
    info = store.verify(args.version)
    payload = json.loads(
        (info.path / "candidates.json").read_text(encoding="utf-8"))
    alt_specs = load_links(args.links)
    recount = recount_payload(payload, alt_specs)
    current = {
        cand["variant_id"]: cand.get("molecule_evidence", {})
        for cand in payload.get("candidates", [])
    }
    alternative = recount["molecule_counts"]
    variant_ids = sorted(set(current) | set(alternative))
    changes = []
    for vid in variant_ids:
        before = current.get(vid, {})
        after = alternative.get(vid, {})
        row = {"variant_id": vid}
        changed = False
        for key in ("alt_molecules", "independent_alt_molecules",
                    "shared_ambiguous_molecules",
                    "linkage_resolved_molecules"):
            old_val = before.get(key, 0)
            new_val = after.get(key, 0)
            if old_val != new_val:
                changed = True
            row[key] = {"current": old_val, "alternative": new_val}
        if changed:
            changes.append(row)
    result = {
        "base_version": args.version,
        "alternative_links": args.links,
        "does_not_modify_base_version": True,
        "support_changes": changes,
        "alternative_units": recount["units"],
        "alternative_molecule_conflicts": recount["molecule_conflicts"],
    }
    text = json.dumps(result, indent=2, ensure_ascii=False)
    print(text)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n", encoding="utf-8")
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
        "link-view": _cmd_link_view,
        "link-select": _cmd_link_select,
        "link-withdraw": _cmd_link_withdraw,
        "link-compare": _cmd_link_compare,
    }
    try:
        return handlers[args.command](args)
    except (StoreError, VerificationError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
