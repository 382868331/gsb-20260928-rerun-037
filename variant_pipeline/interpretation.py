"""Provisional interpretation choices and shared-evidence reconciliation."""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .store import RunStore, VersionInfo, sha256_text

SELECTIONS_DIRNAME = "selections"


@dataclass(frozen=True)
class Selection:
    selection_id: str
    base_version: str
    selected_variant_id: str | None
    status: str  # provisional | withdrawn


def selection_path(store: RunStore, version: str) -> Path:
    safe = version.replace("/", "_")
    return store.root / SELECTIONS_DIRNAME / f"{safe}.json"


def save_selection(
    store: RunStore,
    version: str,
    selected_variant_id: str,
) -> Selection:
    info = store.verify(version)
    candidates = {
        c["variant_id"] for c in json.loads(
            (info.path / "candidates.json").read_text(encoding="utf-8")
        ).get("candidates", [])
    }
    if selected_variant_id not in candidates:
        raise ValueError(f"unknown candidate in {version}: {selected_variant_id}")
    target = selection_path(store, version)
    existing = None
    if target.is_file():
        existing = json.loads(target.read_text(encoding="utf-8"))
        if existing.get("status") == "provisional":
            raise ValueError(
                f"{version} already has provisional selection "
                f"{existing.get('selected_variant_id')}; withdraw it first"
            )
    record = {
        "selection_id": existing.get("selection_id") if existing else sha256_text(
            f"{version}:{selected_variant_id}")[:16],
        "base_version": version,
        "selected_variant_id": selected_variant_id,
        "status": "provisional",
        "supersedes": existing.get("selection_id") if existing else None,
    }
    _write_json(target, record)
    return Selection(
        record["selection_id"], version, selected_variant_id, "provisional")


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
        record["selection_id"], version,
        record.get("selected_variant_id"), "withdrawn")


def load_selection(store: RunStore, version: str) -> dict | None:
    path = selection_path(store, version)
    if not path.is_file():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("status") == "provisional" else None


def selected_view(info: VersionInfo, selected_variant_id: str) -> dict:
    payload = json.loads(
        (info.path / "candidates.json").read_text(encoding="utf-8"))
    selected_group = next(
        (group for group in payload.get("ambiguity_groups", [])
         if selected_variant_id in group.get("candidate_variant_ids", [])),
        None,
    )
    assigned_ids: set[str] = set()
    if selected_group is not None:
        assigned_ids = {
            record["fragment_id"]
            for record in selected_group["shared_fragments"]
            if selected_variant_id in record.get("supports", [])
        }
    candidates = []
    for candidate in payload["candidates"]:
        candidate = json.loads(json.dumps(candidate))
        if selected_group is None:
            candidate["selection_state"] = (
                "provisionally_selected"
                if candidate["variant_id"] == selected_variant_id
                else "unchanged"
            )
            candidates.append(candidate)
            continue
        if candidate["variant_id"] not in selected_group["candidate_variant_ids"]:
            candidates.append(candidate)
            continue
        support = candidate["supporting_evidence"]
        if candidate["variant_id"] == selected_variant_id:
            candidate["selection_state"] = "provisionally_selected"
            support["provisional_independent_alt_fragments"] = (
                support["independent_alt_fragments"] + len(assigned_ids))
            candidate["status_note"] = (
                "shared ambiguous fragments are assigned only in this temporary view")
        else:
            candidate["selection_state"] = "deprioritized_by_provisional_choice"
            retained_shared = sum(
                1
                for record in selected_group["shared_fragments"]
                if candidate["variant_id"] in record.get("supports", [])
                and selected_variant_id not in record.get("supports", [])
            )
            support["provisional_independent_alt_fragments"] = (
                support["independent_alt_fragments"] + retained_shared)
            candidate["status_note"] = (
                "shared ambiguous fragments are assigned to the selected candidate "
                "in this temporary view and add no evidence here")
        support["provisional_assigned_shared_fragment_count"] = len(assigned_ids)
        candidates.append(candidate)
    payload["candidates"] = candidates
    payload["provisional_selection"] = {
        "base_version": info.version,
        "selected_variant_id": selected_variant_id,
        "does_not_modify_base_version": True,
        "ambiguity_group_id": (
            selected_group["group_id"] if selected_group else None),
        "assigned_fragment_ids": sorted(assigned_ids),
    }
    return payload


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".selection-", dir=str(path.parent))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
