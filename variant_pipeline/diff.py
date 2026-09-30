"""Compare two frozen batch results.

The diff reads stored candidates.json payloads only. It never re-analyzes
historical inputs, so candidate status and evidence counts cannot be
rewritten by newer code. Each per-candidate change is attributed to either
an input change (different reads/reference/config snapshot) or an
interpretation-rule change (different analysis_version), using the input
hashes recorded in each batch manifest. Support/opposition evidence is read
from each batch as-is and displayed side by side: records from two batches
are never merged or concatenated.
"""

from __future__ import annotations

from . import store

TRACKED_FIELDS = (
    "alt_fragments",
    "low_quality_alt_fragments",
    "reference_fragments",
    "other_allele_fragments",
    "overlap_conflict_fragments",
    "depth",
    "passing_depth",
    "informative_depth",
    "confidence",
    "status",
)


def _inputs_changed(manifest_a: dict, manifest_b: dict) -> bool:
    return manifest_a.get("input_hashes") != manifest_b.get("input_hashes")


def _rules_changed(manifest_a: dict, manifest_b: dict) -> bool:
    return manifest_a.get("analysis_version") != manifest_b.get("analysis_version")


def _candidate_map(payload: dict) -> dict[str, dict]:
    return {c["variant_id"]: c for c in payload.get("candidates", [])}


def _field_values(cand: dict, field: str):
    if field in {
        "low_quality_alt_fragments",
        "depth",
        "passing_depth",
        "informative_depth",
        "confidence",
        "status",
    }:
        return cand.get(field)
    supporting = cand.get("supporting_evidence", {})
    opposing = cand.get("opposing_evidence", {})
    if field == "alt_fragments":
        return supporting.get("alt_fragments")
    if field == "reference_fragments":
        return opposing.get("reference_fragments")
    if field == "other_allele_fragments":
        return opposing.get("other_allele_fragments")
    if field == "overlap_conflict_fragments":
        return opposing.get("overlap_conflict_fragments")
    return None


def _attribution(inputs_changed: bool, rules_changed: bool) -> list[str]:
    causes = []
    if inputs_changed:
        causes.append("input_change")
    if rules_changed:
        causes.append("interpretation_rule_change")
    return causes or ["unchanged_inputs_and_rules"]


def diff_batches(root, batch_a: str, batch_b: str) -> dict:
    bid_a = store.resolve_batch_id(root, batch_a)
    bid_b = store.resolve_batch_id(root, batch_b)
    manifest_a = store.load_manifest(root, bid_a)
    manifest_b = store.load_manifest(root, bid_b)
    payload_a = store.load_payload(root, bid_a)
    payload_b = store.load_payload(root, bid_b)

    inputs_changed = _inputs_changed(manifest_a, manifest_b)
    rules_changed = _rules_changed(manifest_a, manifest_b)
    causes = _attribution(inputs_changed, rules_changed)

    map_a = _candidate_map(payload_a)
    map_b = _candidate_map(payload_b)
    added, removed, changed = [], [], []

    for vid in sorted(set(map_a) | set(map_b)):
        in_a, in_b = vid in map_a, vid in map_b
        if in_a and not in_b:
            removed.append(
                {
                    "variant_id": vid,
                    "change": "disappeared",
                    "attribution": causes,
                    "from": _candidate_snapshot(map_a[vid]),
                    "to": None,
                }
            )
        elif in_b and not in_a:
            added.append(
                {
                    "variant_id": vid,
                    "change": "appeared",
                    "attribution": causes,
                    "from": None,
                    "to": _candidate_snapshot(map_b[vid]),
                }
            )
        else:
            field_changes = []
            for field in TRACKED_FIELDS:
                old = _field_values(map_a[vid], field)
                new = _field_values(map_b[vid], field)
                if old != new:
                    field_changes.append({"field": field, "from": old, "to": new})
            if field_changes:
                changed.append(
                    {
                        "variant_id": vid,
                        "change": "modified",
                        "attribution": causes,
                        "field_changes": field_changes,
                        "from": _candidate_snapshot(map_a[vid]),
                        "to": _candidate_snapshot(map_b[vid]),
                    }
                )

    return {
        "schema": "variant-pipeline-diff/1",
        "from_batch": {
            "batch_id": bid_a,
            "analysis_version": manifest_a["analysis_version"],
            "input_hashes": manifest_a["input_hashes"],
        },
        "to_batch": {
            "batch_id": bid_b,
            "analysis_version": manifest_b["analysis_version"],
            "input_hashes": manifest_b["input_hashes"],
        },
        "inputs_changed": inputs_changed,
        "interpretation_rules_changed": rules_changed,
        "added": added,
        "removed": removed,
        "changed": changed,
    }


def _candidate_snapshot(cand: dict) -> dict:
    supporting = cand.get("supporting_evidence", {})
    opposing = cand.get("opposing_evidence", {})
    return {
        "variant_id": cand["variant_id"],
        "status": cand.get("status"),
        "confidence": cand.get("confidence"),
        "alt_fragments": supporting.get("alt_fragments"),
        "low_quality_alt_fragments": supporting.get(
            "low_quality_alt_fragments"
        ),
        "reference_fragments": opposing.get("reference_fragments"),
        "other_allele_fragments": opposing.get("other_allele_fragments"),
        "overlap_conflict_fragments": opposing.get(
            "overlap_conflict_fragments"
        ),
        "depth": cand.get("depth"),
        "flags": list(cand.get("flags", [])),
    }


def _cause_label(causes: list[str]) -> str:
    labels = {
        "input_change": "输入变化",
        "interpretation_rule_change": "方向解释规则变化",
        "unchanged_inputs_and_rules": "输入与规则均未变(重复批次)",
    }
    return " + ".join(labels[c] for c in causes)


def render_diff_text(result: dict) -> str:
    lines: list[str] = []
    lines.append("批次间变异差异比较 (research use only / 仅供研究流程)")
    lines.append("=" * 68)
    lines.append(f"旧批次: {result['from_batch']['batch_id']}")
    lines.append(f"新批次: {result['to_batch']['batch_id']}")
    lines.append(
        f"输入变化: {'是' if result['inputs_changed'] else '否'}  "
        f"方向解释规则变化: "
        f"{'是' if result['interpretation_rules_changed'] else '否'}"
    )
    lines.append("")

    def describe(entry: dict) -> None:
        cause = _cause_label(entry["attribution"])
        lines.append(f"- [{entry['change']}] {entry['variant_id']}  归因: {cause}")
        for snap_key, label in (("from", "旧"), ("to", "新")):
            snap = entry.get(snap_key)
            if snap is None:
                continue
            lines.append(
                f"    {label}: 状态={snap['status']} 可信度={snap['confidence']} "
                f"ALT={snap['alt_fragments']} REF={snap['reference_fragments']} "
                f"其他等位={snap['other_allele_fragments']} "
                f"重叠冲突={snap['overlap_conflict_fragments']} 深度={snap['depth']}"
            )
        for change in entry.get("field_changes", []):
            lines.append(
                f"    {change['field']}: {change['from']} -> {change['to']}"
            )

    lines.append(f"新增候选 ({len(result['added'])})")
    lines.append("-" * 68)
    for entry in result["added"]:
        describe(entry)
    if not result["added"]:
        lines.append("无")
    lines.append("")
    lines.append(f"消失候选 ({len(result['removed'])})")
    lines.append("-" * 68)
    for entry in result["removed"]:
        describe(entry)
    if not result["removed"]:
        lines.append("无")
    lines.append("")
    lines.append(f"支持数量/状态变化候选 ({len(result['changed'])})")
    lines.append("-" * 68)
    for entry in result["changed"]:
        describe(entry)
    if not result["changed"]:
        lines.append("无")
    lines.append("")
    lines.append(
        "注: 两批次证据各自独立保留(含共享片段与软截断的原始记录), 未做拼接。"
    )
    return "\n".join(lines) + "\n"
