"""Compare two published analysis versions candidate by candidate.

The comparison attributes every difference to what actually changed
between the two runs: the input data (reads / reference / sample) or the
interpretation rules (analysis version / parameters). Evidence records are
never merged across versions: for a candidate present in both, the report
shows each version's own support counts side by side, so shared fragments
and soft-clipped evidence keep their original, per-version records and
are never spliced into a combined tally.
"""

from __future__ import annotations

import json

from .store import RunStore, VersionInfo

SUPPORT_FIELDS = (
    "alt_fragments",
    "forward",
    "reverse",
    "low_quality_alt_fragments",
)

OPPOSE_FIELDS = (
    "reference_fragments",
    "other_allele_fragments",
    "overlap_conflict_fragments",
)


def _load_payload(info: VersionInfo) -> dict:
    return json.loads(
        (info.path / "candidates.json").read_text(encoding="utf-8"))


def _candidate_map(payload: dict) -> dict:
    return {c["variant_id"]: c for c in payload.get("candidates", [])}


def attribute_cause(from_info: VersionInfo, to_info: VersionInfo) -> dict:
    """Classify what changed between the two runs (input vs rules)."""
    from_inputs = from_info.manifest.get("inputs", {})
    to_inputs = to_info.manifest.get("inputs", {})
    input_changed = any(
        from_inputs.get(key) != to_inputs.get(key)
        for key in ("reads_sha256", "fasta_sha256")
    )
    from_params = from_inputs.get("config", {}).get("parameters", {})
    to_params = to_inputs.get("config", {}).get("parameters", {})
    rule_changed = (
        from_info.manifest.get("analysis_version")
        != to_info.manifest.get("analysis_version")
        or from_params != to_params
    )
    if input_changed and rule_changed:
        cause = "input_and_rule_change"
    elif input_changed:
        cause = "input_change"
    elif rule_changed:
        cause = "rule_change"
    else:
        cause = "no_input_or_rule_change"
    return {
        "cause": cause,
        "input_changed": input_changed,
        "rule_changed": rule_changed,
        "from_analysis_version": from_info.manifest.get("analysis_version"),
        "to_analysis_version": to_info.manifest.get("analysis_version"),
        "from_parameters": from_params,
        "to_parameters": to_params,
    }


def _evidence_snapshot(candidate: dict) -> dict:
    """Per-version evidence snapshot; never combined across versions."""
    support = candidate.get("supporting_evidence", {})
    oppose = candidate.get("opposing_evidence", {})
    return {
        "status": candidate.get("status"),
        "confidence": candidate.get("confidence"),
        "support": {k: support.get(k) for k in SUPPORT_FIELDS},
        "oppose": {k: oppose.get(k) for k in OPPOSE_FIELDS},
    }


def _changed_fields(from_cand: dict, to_cand: dict) -> dict:
    from_snap = _evidence_snapshot(from_cand)
    to_snap = _evidence_snapshot(to_cand)
    changes = {}
    for key in ("status", "confidence"):
        if from_snap[key] != to_snap[key]:
            changes[key] = {"from": from_snap[key], "to": to_snap[key]}
    for section in ("support", "oppose"):
        for key in from_snap[section]:
            from_val = from_snap[section][key]
            to_val = to_snap[section][key]
            if from_val != to_val:
                changes[f"{section}.{key}"] = {"from": from_val, "to": to_val}
    return changes


def compare_versions(store: RunStore, from_version: str,
                     to_version: str) -> dict:
    """Diff two published versions; both are verified before reading."""
    from_info = store.verify(from_version)
    to_info = store.verify(to_version)
    from_cands = _candidate_map(_load_payload(from_info))
    to_cands = _candidate_map(_load_payload(to_info))
    attribution = attribute_cause(from_info, to_info)

    added, removed, changed = [], [], []
    for variant_id in sorted(to_cands):
        if variant_id not in from_cands:
            added.append({
                "variant_id": variant_id,
                "cause": attribution["cause"],
                "to": _evidence_snapshot(to_cands[variant_id]),
            })
    for variant_id in sorted(from_cands):
        if variant_id not in to_cands:
            removed.append({
                "variant_id": variant_id,
                "cause": attribution["cause"],
                "from": _evidence_snapshot(from_cands[variant_id]),
            })
    for variant_id in sorted(from_cands.keys() & to_cands.keys()):
        fields = _changed_fields(from_cands[variant_id], to_cands[variant_id])
        if fields:
            changed.append({
                "variant_id": variant_id,
                "cause": attribution["cause"],
                "changes": fields,
                # Original per-version evidence, kept side by side and
                # never summed or merged across versions.
                "from": _evidence_snapshot(from_cands[variant_id]),
                "to": _evidence_snapshot(to_cands[variant_id]),
            })
    return {
        "from_version": from_version,
        "to_version": to_version,
        "attribution": attribution,
        "added": added,
        "removed": removed,
        "changed": changed,
        "summary": {
            "added": len(added),
            "removed": len(removed),
            "changed": len(changed),
            "shared_unchanged": len(from_cands.keys() & to_cands.keys())
            - len(changed),
        },
    }


def render_diff_text(diff: dict) -> str:
    attr = diff["attribution"]
    cause_labels = {
        "input_change": "输入变化",
        "rule_change": "方向/规则解释变化",
        "input_and_rule_change": "输入与规则同时变化",
        "no_input_or_rule_change": "输入与规则均未变化",
    }
    cause = cause_labels.get(attr["cause"], attr["cause"])
    lines = [
        f"版本比较: {diff['from_version']} -> {diff['to_version']}",
        "=" * 68,
        f"归因: {cause} "
        f"(输入变化={attr['input_changed']}, 规则变化={attr['rule_changed']})",
        f"规则版本: {attr['from_analysis_version']} -> "
        f"{attr['to_analysis_version']}",
        "说明: 共享候选的支持/反对证据按各自版本分别列出，"
        "不跨版本合并或拼接。",
        "",
    ]
    s = diff["summary"]
    lines.append(
        f"汇总: 新增={s['added']} 消失={s['removed']} "
        f"支持变化={s['changed']} 不变={s['shared_unchanged']}"
    )
    lines.append("")
    for entry in diff["added"]:
        lines.append(f"[新增] {entry['variant_id']} (归因: {cause})")
        lines.append(f"  新版本证据: {json.dumps(entry['to'], ensure_ascii=False)}")
    for entry in diff["removed"]:
        lines.append(f"[消失] {entry['variant_id']} (归因: {cause})")
        lines.append(f"  旧版本证据: {json.dumps(entry['from'], ensure_ascii=False)}")
    for entry in diff["changed"]:
        lines.append(f"[变化] {entry['variant_id']} (归因: {cause})")
        for field, vals in entry["changes"].items():
            lines.append(f"  {field}: {vals['from']} -> {vals['to']}")
    return "\n".join(lines) + "\n"

