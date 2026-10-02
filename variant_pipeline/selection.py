"""Tentative interpretation choices kept outside immutable run outputs."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

from .store import StoreError

SELECTION_FILENAME = "selections.json"


def selection_path(store) -> Path:
    return Path(store.root) / SELECTION_FILENAME


def load_selections(store) -> dict:
    path = selection_path(store)
    if not path.is_file():
        return {"schema": "interpretation-selections-v1", "versions": {}}
    return json.loads(path.read_text(encoding="utf-8"))


def save_selections(store, data: dict) -> None:
    store.root.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".selections-", dir=str(store.root))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, selection_path(store))


def load_payload(store, version: str) -> dict:
    info = store.verify(version)
    return json.loads((info.path / "candidates.json").read_text(encoding="utf-8"))


def ambiguity_groups(payload: dict) -> list[dict]:
    known = {c["variant_id"] for c in payload.get("candidates", [])}
    groups = []
    for record in payload.get("evidence_ledger", {}).get("fragment_records", []):
        for raw_ids in record.get("ambiguity_groups", []):
            ids = [v for v in raw_ids if v in known]
            if len(ids) > 1:
                groups.append({
                    "fragment_uid": record["fragment_uid"],
                    "variant_ids": sorted(ids),
                    "group_id": "|".join(sorted(ids)),
                })
    return groups


def _version_decisions(data: dict, version: str) -> dict:
    return data.setdefault("versions", {}).setdefault(version, {})


def adopt(store, version: str, variant_id: str) -> dict:
    payload = load_payload(store, version)
    groups = ambiguity_groups(payload)
    containing = [g for g in groups if variant_id in g["variant_ids"]]
    if not containing:
        raise StoreError(
            f"{variant_id} is not part of a shared ambiguity group in {version}")
    data = load_selections(store)
    decisions = _version_decisions(data, version)
    adopted = decisions.setdefault("adopted", {})
    for group in containing:
        current = adopted.get(group["group_id"])
        if current == variant_id:
            continue
        if current is not None:
            raise StoreError(
                f"manual choice already exists for this ambiguity: {current}")
        adopted[group["group_id"]] = variant_id
    save_selections(store, data)
    return data["versions"][version]


def withdraw(store, version: str, variant_id: str) -> dict:
    data = load_selections(store)
    decisions = data.get("versions", {}).get(version, {})
    adopted = decisions.get("adopted", {})
    decisions["adopted"] = {
        group_id: chosen for group_id, chosen in adopted.items()
        if chosen != variant_id
    }
    save_selections(store, data)
    return decisions


def selected_variant_ids(store, version: str) -> set[str]:
    data = load_selections(store)
    decisions = data.get("versions", {}).get(version, {})
    return set(decisions.get("adopted", {}).values())


def project(payload: dict, selected: set[str]) -> dict:
    groups = ambiguity_groups(payload)
    candidates = {c["variant_id"]: c for c in payload.get("candidates", [])}
    assigned = {vid: [] for vid in candidates}
    withheld = {vid: [] for vid in candidates}
    for group in groups:
        chosen = [vid for vid in group["variant_ids"] if vid in selected]
        if len(chosen) == 1:
            chosen_id = chosen[0]
            assigned[chosen_id].append(group["fragment_uid"])
            for vid in group["variant_ids"]:
                if vid != chosen_id:
                    withheld[vid].append(group["fragment_uid"])
        elif chosen:
            # Multiple chosen members remain mutually incompatible; no
            # shared fragment can count as independent evidence.
            for vid in group["variant_ids"]:
                withheld[vid].append(group["fragment_uid"])

    projected_candidates = []
    for vid, candidate in sorted(candidates.items()):
        support = candidate["supporting_evidence"]
        shared = support.get("shared_ambiguous_alt_fragments", 0)
        lost = len(withheld[vid])
        projected_candidates.append({
            "variant_id": vid,
            "baseline_alt_fragments": support.get("alt_fragments"),
            "independent_alt_fragments": support.get(
                "independent_alt_fragments", 0),
            "projected_alt_fragments": max(
                0,
                support.get("independent_alt_fragments", 0)
                + shared
                - lost,
            ),
            "assigned_shared_fragments": sorted(assigned[vid]),
            "withheld_shared_fragments": sorted(withheld[vid]),
        })
    return {
        "schema": "tentative-interpretation-projection-v1",
        "selected_variant_ids": sorted(selected),
        "candidates": projected_candidates,
    }


def render_projection(projection: dict) -> str:
    lines = [
        "临时解释投影（不修改原始分析版本）",
        "=" * 68,
        "已选: " + ", ".join(projection["selected_variant_ids"] or ["无"]),
        "",
    ]
    for cand in projection["candidates"]:
        lines.append(
            f"{cand['variant_id']}: "
            f"基线ALT={cand['baseline_alt_fragments']} -> "
            f"投影ALT={cand['projected_alt_fragments']} "
            f"(独立={cand['independent_alt_fragments']}, "
            f"分得共享={len(cand['assigned_shared_fragments'])}, "
            f"移出共享={len(cand['withheld_shared_fragments'])})"
        )
    return "\n".join(lines) + "\n"
