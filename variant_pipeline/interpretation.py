"""Provisional interpretation choices and evidence reconciliation."""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .store import RunStore, VersionInfo, sha256_text

SELECTIONS_DIRNAME = "selections"
_HISTORY_FIELDS = (
    "selection_id", "selected_variant_id", "evidence_unit_id", "relation_id",
    "relation_version", "rationale", "excluded_variant_ids", "adopted_basis",
    "conflicting_variant_ids", "conflicting_combinations",
    "created_at", "superseded_at",
)


@dataclass(frozen=True)
class Selection:
    selection_id: str
    base_version: str
    selected_variant_id: str | None
    status: str


def selection_path(store: RunStore, version: str) -> Path:
    safe = version.replace("/", "_")
    return store.root / SELECTIONS_DIRNAME / f"{safe}.json"


def _payload(info: VersionInfo) -> dict:
    return json.loads((info.path / "candidates.json").read_text(encoding="utf-8"))


def save_selection(
    store: RunStore,
    version: str,
    selected_variant_id: str,
    evidence_unit_id: str | None = None,
    relation_id: str | None = None,
    rationale: str | None = None,
) -> Selection:
    info = store.verify(version)
    payload = _payload(info)
    candidates = {c["variant_id"] for c in payload.get("candidates", [])}
    if selected_variant_id not in candidates:
        raise ValueError(f"unknown candidate in {version}: {selected_variant_id}")

    target = selection_path(store, version)
    existing = json.loads(target.read_text(encoding="utf-8")) if target.is_file() else None
    history = list(existing.get("history", [])) if existing else []
    if existing and existing.get("status") == "provisional":
        archived = {key: existing.get(key) for key in _HISTORY_FIELDS}
        history.append(archived)

    unit = None
    if evidence_unit_id:
        unit = next((u for u in payload.get("evidence_units", [])
                     if u.get("association_id") == evidence_unit_id), None)
        if unit is None:
            raise ValueError(f"unknown evidence unit in {version}: {evidence_unit_id}")
        unit_candidates = set(unit.get("linked_candidates", [])) | set(
            unit.get("conflicting_candidates", []))
        if selected_variant_id not in unit_candidates:
            raise ValueError(
                f"selected variant {selected_variant_id} is not observed by "
                f"evidence unit {evidence_unit_id}")

    selected_group = None
    excluded: list[str] = []
    conflicting: list[str] = []
    conflict_basis: list[dict] = []
    adopted_basis = "ambiguity_group"
    if unit is not None:
        single = {
            candidate
            for values in unit.get("single_fragment_candidates", {}).values()
            for candidate in values
        }
        conflicting = list(unit.get("conflicting_candidates", []))
        conflict_basis = list(unit.get("conflicting_combinations", []))
        excluded = sorted(
            single
            - set(unit.get("linked_candidates", []))
            - set(conflicting))
        adopted_basis = "evidence_unit_relation"
    else:
        selected_group = next(
            (group for group in payload.get("ambiguity_groups", [])
             if selected_variant_id in group.get("candidate_variant_ids", [])),
            None,
        )
        if selected_group is not None:
            excluded = sorted(
                set(selected_group.get("candidate_variant_ids", []))
                - {selected_variant_id}
            )

    record = {
        "selection_id": sha256_text(
            f"{version}:{selected_variant_id}:{relation_id}:{evidence_unit_id}:{len(history)}"
        )[:16],
        "base_version": version,
        "selected_variant_id": selected_variant_id,
        "status": "provisional",
        "evidence_unit_id": evidence_unit_id,
        "relation_id": relation_id or (unit or {}).get("relation_id"),
        "relation_version": (unit or {}).get("relation_version"),
        "rationale": rationale,
        "excluded_variant_ids": excluded,
        "conflicting_variant_ids": conflicting,
        "conflicting_combinations": conflict_basis,
        "adopted_basis": adopted_basis,
        "supersedes": existing.get("selection_id") if existing else None,
        "history": history,
    }
    _write_json(target, record)
    return Selection(record["selection_id"], version, selected_variant_id, "provisional")


def withdraw_selection(store: RunStore, version: str) -> Selection:
    target = selection_path(store, version)
    if not target.is_file():
        raise ValueError(f"no provisional selection for {version}")
    record = json.loads(target.read_text(encoding="utf-8"))
    if record.get("status") == "withdrawn":
        raise ValueError(f"selection for {version} is already withdrawn")
    record["status"] = "withdrawn"
    _write_json(target, record)
    return Selection(
        record["selection_id"], version, record.get("selected_variant_id"), "withdrawn")


def load_selection(store: RunStore, version: str) -> dict | None:
    path = selection_path(store, version)
    if not path.is_file():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("status") == "provisional" else None


def selected_view(info: VersionInfo, selected_variant_id: str,
                  decision: dict | None = None) -> dict:
    payload = _payload(info)
    decision = decision or {}
    selected_group = next(
        (group for group in payload.get("ambiguity_groups", [])
         if selected_variant_id in group.get("candidate_variant_ids", [])),
        None,
    )
    assigned_ids: set[str] = set()
    if selected_group is not None:
        assigned_ids = {
            record["fragment_id"]
            for record in selected_group.get("shared_fragments", [])
            if selected_variant_id in record.get("supports", [])
        }

    selected_unit = None
    unit_ids: set[str] = set()
    if decision.get("evidence_unit_id"):
        selected_unit = next(
            (unit for unit in payload.get("evidence_units", [])
             if unit.get("association_id") == decision.get("evidence_unit_id")),
            None,
        )
        if selected_unit is not None:
            unit_ids = {
                variant_id
                for values in selected_unit.get("single_fragment_candidates", {}).values()
                for variant_id in values
            }

    candidates = []
    group_ids = set(selected_group.get("candidate_variant_ids", [])) if selected_group else set()
    excluded_by_decision = set(decision.get("excluded_variant_ids", []))
    for original in payload["candidates"]:
        candidate = json.loads(json.dumps(original))
        variant_id = candidate["variant_id"]

        if selected_unit is not None and variant_id in unit_ids:
            if variant_id == selected_variant_id:
                candidate["selection_state"] = "provisionally_selected"
                candidate["status_note"] = (
                    "selected linked evidence-unit interpretation; "
                    "the unit remains one independent evidence source")
            elif variant_id in decision.get("conflicting_variant_ids", []):
                candidate["selection_state"] = "intra_molecule_conflict"
                candidate["status_note"] = (
                    "retained allele contradicts the selected allele within "
                    "the same molecule; not excluded by import order")
            elif variant_id in excluded_by_decision:
                candidate["selection_state"] = "excluded_by_provisional_choice"
                candidate["status_note"] = "excluded only under the recorded relation version"
            else:
                candidate["selection_state"] = "linked_support_unchanged"
            candidates.append(candidate)
            continue

        if variant_id not in group_ids:
            candidate["selection_state"] = (
                "provisionally_selected" if variant_id == selected_variant_id else "unchanged")
            candidates.append(candidate)
            continue

        support = candidate["supporting_evidence"]
        if variant_id == selected_variant_id:
            candidate["selection_state"] = "provisionally_selected"
            support["provisional_independent_alt_fragments"] = (
                support["independent_alt_fragments"] + len(assigned_ids))
            candidate["status_note"] = "ambiguous fragment assigned only in this temporary view"
        else:
            candidate["selection_state"] = "deprioritized_by_provisional_choice"
            retained_shared = sum(
                1
                for record in selected_group["shared_fragments"]
                if variant_id in record.get("supports", [])
                and selected_variant_id not in record.get("supports", []))
            support["provisional_independent_alt_fragments"] = (
                support["independent_alt_fragments"] + retained_shared)
            candidate["status_note"] = "assigned ambiguous fragment adds no evidence here"
        support["provisional_assigned_shared_fragment_count"] = len(assigned_ids)
        candidates.append(candidate)

    excluded = excluded_by_decision if decision else group_ids - {selected_variant_id}
    payload["candidates"] = candidates
    payload["provisional_selection"] = {
        "base_version": info.version,
        "selected_variant_id": selected_variant_id,
        "does_not_modify_base_version": True,
        "ambiguity_group_id": selected_group["group_id"] if selected_group else None,
        "assigned_fragment_ids": sorted(assigned_ids),
        "excluded_variant_ids": sorted(excluded),
        "conflicting_variant_ids": decision.get("conflicting_variant_ids", []),
        "conflicting_combinations": decision.get("conflicting_combinations", []),
        "evidence_unit_id": decision.get("evidence_unit_id"),
        "relation_id": decision.get("relation_id"),
        "relation_version": decision.get("relation_version"),
        "adopted_basis": decision.get("adopted_basis", "ambiguity_group"),
        "rationale": decision.get("rationale"),
    }
    return payload


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".selection-", dir=str(path.parent))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
