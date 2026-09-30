"""Stable, annotation-ready output writers.

Each candidate carries the canonical variant identity (``variant_id`` plus
chrom/pos/ref/alt) so a later annotation stage can join on the same key.
"""

from __future__ import annotations

import csv
import io
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
        "alt_fragments": ev.alt_fragments,
        "ref_fragments": ev.ref_fragments,
        "other_alt_fragments": ev.other_alt_fragments,
        "conflicting_overlap_fragments": ev.conflicting_overlap_fragments,
        "depth": ev.depth,
        "passing_depth": ev.passing_depth,
        "informative_depth": ev.informative_depth,
        "supporting_evidence": ev.supporting_summary(),
        "opposing_evidence": ev.opposing_summary(),
        "coverage": {
            "fragment_depth": ev.depth,
            "passing_depth": ev.passing_depth,
            "informative_depth": ev.informative_depth,
        },
    }


def build_payload(cfg: RunConfig, result: CallResult, batch: dict | None = None) -> dict:
    payload = {
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
    }
    if batch is not None:
        payload["batch"] = batch
    return payload


def write_json(payload: dict, path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def write_candidates_tsv(result: CallResult, path: str | Path) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        handle.write(candidates_tsv_text(result))


def candidates_tsv_text(result: CallResult) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=TSV_COLUMNS, delimiter="\t")
    writer.writeheader()
    for ev in result.candidates:
        row = {col: getattr(ev, col, "") for col in TSV_COLUMNS if col != "flags"}
        row["flags"] = ";".join(ev.flags)
        writer.writerow(row)
    return buffer.getvalue()


def write_coverage_tsv(result: CallResult, path: str | Path) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        handle.write(coverage_tsv_text(result))


def coverage_tsv_text(result: CallResult) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t")
    writer.writerow(["chrom", "pos", "fragment_depth", "passing_depth"])
    for c in result.coverage:
        writer.writerow([c.chrom, c.pos, c.fragment_depth, c.passing_depth])
    return buffer.getvalue()


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
        lines.append(
            f"[{cand['status']}] {cand['variant_id']} ({cand['kind']}) "
            f"可信度={cand['confidence']}"
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
