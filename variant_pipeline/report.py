"""Stable, annotation-ready output writers.

Each candidate carries the canonical variant identity (``variant_id`` plus
chrom/pos/ref/alt) so a later annotation stage can join on the same key.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from . import ANALYSIS_VERSION, RESEARCH_DISCLAIMER
from .caller import CallResult
from .config import RunConfig


TSV_COLUMNS = [
    "variant_id", "chrom", "pos", "ref", "alt", "kind", "status",
    "confidence", "alt_fragments", "ref_fragments", "other_alt_fragments",
    "alt_forward", "alt_reverse", "strand_bias_p",
    "low_quality_alt_fragments", "conflicting_overlap_fragments",
    "independent_alt_fragments", "shared_ambiguous_alt_fragments",
    "alt_molecules", "independent_alt_molecules",
    "linkage_resolved_molecules", "shared_ambiguous_molecules",
    "depth", "passing_depth", "informative_depth", "alt_fraction",
    "alt_mean_base_quality", "flags",
]


def candidate_dict(ev) -> dict:
    return {
        "variant_id": ev.variant_id,
        "chrom": ev.chrom,
        "pos": ev.pos,
        "ref": ev.ref,
        "alt": ev.alt,
        "kind": ev.kind,
        "status": ev.status,
        "confidence": ev.confidence,
        "score_components": ev.components,
        "flags": ev.flags,
        "supporting_evidence": ev.supporting_summary(),
        "evidence_ambiguity": {
            "independent_alt_fragments": ev.independent_alt_fragments,
            "shared_ambiguous_alt_fragments": ev.shared_ambiguous_alt_fragments,
            "shared_fragments": [
                {"fragment_id": fid, "also_supports": vids}
                for fid, vids in sorted(ev.shared_with.items())
            ],
        },
        "molecule_evidence": ev.molecule_summary(),
        "opposing_evidence": ev.opposing_summary(),
        "coverage": {
            "fragment_depth": ev.depth,
            "passing_depth": ev.passing_depth,
            "informative_depth": ev.informative_depth,
        },
    }


def build_payload(cfg: RunConfig, result: CallResult) -> dict:
    return {
        "disclaimer": RESEARCH_DISCLAIMER,
        "analysis_version": ANALYSIS_VERSION,
        "run": {
            "sample": {
                "sample_id": cfg.sample.sample_id,
                "read_group": cfg.sample.read_group,
                "description": cfg.sample.description,
            },
            "reference": {
                "name": cfg.reference.name,
                "assembly": cfg.reference.assembly,
                "fasta_path": cfg.reference.fasta_path,
            },
            "parameters": cfg.parameters.__dict__,
            "reads_path": cfg.reads_path,
        },
        "candidates": [candidate_dict(ev) for ev in result.candidates],
        "low_coverage_regions": [
            {"chrom": c, "start": s, "end": e, "min_fragment_depth": d}
            for c, s, e, d in result.low_coverage
        ],
        "ambiguity_groups": result.ambiguity_groups,
        "fragment_sources": result.fragment_sources,
        "linkage_units": result.linkage_units,
        "molecule_conflicts": result.molecule_conflicts,
        "fragment_linkage_detail": result.fragment_linkage_detail,
    }


def write_json(payload: dict, path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def write_candidates_tsv(result: CallResult, path: str | Path) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for ev in result.candidates:
            row = {col: getattr(ev, col, "") for col in TSV_COLUMNS if col != "flags"}
            row["flags"] = ";".join(ev.flags)
            writer.writerow(row)


def write_coverage_tsv(result: CallResult, path: str | Path) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chrom", "pos", "fragment_depth", "passing_depth"])
        for c in result.coverage:
            writer.writerow([c.chrom, c.pos, c.fragment_depth, c.passing_depth])


def render_text_report(payload: dict) -> str:
    lines: list[str] = []
    run = payload["run"]
    lines.append("候选变异证据汇总 (research use only / 仅供研究流程)")
    lines.append("=" * 68)
    lines.append(f"分析规则版本: {payload['analysis_version']}")
    lines.append(f"样本: {run['sample']['sample_id']}  参考版本: {run['reference']['assembly']}")
    lines.append("提示: " + payload["disclaimer"])
    lines.append("")
    for cand in payload["candidates"]:
        s = cand["supporting_evidence"]
        o = cand["opposing_evidence"]
        cov = cand["coverage"]
        amb = cand.get("evidence_ambiguity", {})
        lines.append(
            f"[{cand['status']}] {cand['variant_id']} ({cand['kind']}) "
            f"可信度={cand['confidence']}"
        )
        lines.append(
            f"  证据归属: 独立支持={amb.get('independent_alt_fragments', 0)} "
            f"共享歧义={amb.get('shared_ambiguous_alt_fragments', 0)}; "
            "共享片段不得重复计为新增独立证据"
        )
        mol = cand.get("molecule_evidence", {})
        lines.append(
            f"  分子单元: 支持分子={mol.get('alt_molecules', 0)} "
            f"独立={mol.get('independent_alt_molecules', 0)} "
            f"关联解决={mol.get('linkage_resolved_molecules', 0)} "
            f"仍歧义={mol.get('shared_ambiguous_molecules', 0)}; "
            "同一原始分子的多个片段只算一个证据单元"
        )
        lines.append(
            f"  支持: ALT片段={s['alt_fragments']} "
            f"(正向={s['forward']}, 反向={s['reverse']}) "
            f"ALT比例={s['alt_fraction_of_informative']} "
            f"平均碱基质量={s['mean_base_quality']} "
            f"低质量ALT片段={s['low_quality_alt_fragments']}"
        )
        lines.append(
            f"  反对: REF片段={o['reference_fragments']} "
            f"其他等位片段={o['other_allele_fragments']} "
            f"重叠矛盾片段={o['overlap_conflict_fragments']} "
            f"方向偏差Fisher p={o['strand_bias_fisher_p']}"
        )
        lines.append(
            f"  覆盖: 片段深度={cov['fragment_depth']} "
            f"有效深度={cov['passing_depth']} "
            f"可判定深度={cov['informative_depth']}"
        )
        if cand["flags"]:
            lines.append("  标记: " + ", ".join(cand["flags"]))
        lines.append("")
    low = payload["low_coverage_regions"]
    lines.append("共享歧义片段 (独立支持与共享证据分列)")
    lines.append("-" * 68)
    groups = payload.get("ambiguity_groups", [])
    if not groups:
        lines.append("无")
    else:
        for group in groups:
            lines.append(
                f"  {group['group_id']}: "
                f"候选={', '.join(group['candidate_variant_ids'])}; "
                f"共享片段数={group['shared_fragment_count']} 状态={group['resolution']}"
            )
            for record in group["shared_fragments"]:
                source = next(
                    (s for s in payload.get("fragment_sources", [])
                     if s["fragment_id"] == record["fragment_id"]),
                    {},
                )
                lines.append(
                    f"    {record['fragment_id']} "
                    f"(来源={source.get('source_name', '-')}, qname={source.get('qname', '-')}) "
                    f"-> {', '.join(record['supports'])}"
                )
    lines.append("")
    lines.append("关联片段证据单元 (同一原始分子只计一次)")
    lines.append("-" * 68)
    units = payload.get("linkage_units", [])
    declared = [u for u in units if u.get("link_id")]
    if not declared:
        lines.append("无关联声明；每个片段各自为一个证据单元")
    else:
        for unit in declared:
            members = ", ".join(m["fragment_id"] for m in unit["members"])
            lines.append(
                f"  {unit['unit_id']} (link={unit['link_id']}, "
                f"分子={unit['molecule_id']}) 状态={unit['state']}"
            )
            if members:
                lines.append(f"    成员: {members}")
            if unit["missing_members"]:
                lines.append(
                    f"    缺失成员: {', '.join(unit['missing_members'])}"
                    " (保留待配对，不静默丢弃)"
                )
            for reason in unit["reasons"]:
                lines.append(
                    f"    原因[{reason['kind']}]: {reason['detail']}"
                )
    conflicts = payload.get("molecule_conflicts", [])
    if conflicts:
        lines.append("  分子级冲突 (组合关系排除的解释):")
        for conflict in conflicts:
            lines.append(
                f"    {conflict['unit_id']}: "
                f"候选={', '.join(conflict['candidate_variant_ids'])} "
                f"原因={conflict['reason']}"
            )
    lines.append("")
    lines.append("覆盖不足区间 (单独说明)")
    lines.append("-" * 68)
    if not low:
        lines.append("无")
    else:
        for region in low:
            lines.append(
                f"  {region['chrom']}:{region['start']}-{region['end']} "
                f"(最低片段深度={region['min_fragment_depth']})"
            )
    return "\n".join(lines) + "\n"
