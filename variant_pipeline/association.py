"""Reconcile fragments that are asserted to come from one original molecule.

The association layer never turns two alignments into two independent votes.
It records the single-fragment candidate sets, tests joint placements under an
explicit relation, and keeps late, partial, revised or inconsistent endpoints
visible instead of silently dropping an inconvenient end.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Iterable


PENDING = "pending_endpoint"
MATCH_INSUFFICIENT = "match_insufficient"
RELATION_INCONSISTENT = "relation_inconsistent"
RESOLVED = "resolved"
REVISED = "revised"


@dataclass(frozen=True)
class CandidatePlacement:
    variant_id: str
    chrom: str
    pos: int
    end: int
    strand: str = "F"


@dataclass(frozen=True)
class AssociationRelation:
    association_id: str
    molecule_id: str
    a: str
    b: str
    source_id_a: str
    source_id_b: str
    orientation: str = "any"
    min_distance: int | None = None
    max_distance: int | None = None
    relation_id: str | None = None
    relation_version: str = "1"
    status: str = "active"
    reason: str | None = None


@dataclass
class EvidenceUnitResult:
    association_id: str
    molecule_id: str
    relation_id: str | None
    relation_version: str
    status: str
    fragment_ids: list[str]
    present_fragment_ids: list[str]
    missing_fragment_ids: list[str]
    duplicate_fragment_ids: list[str]
    single_fragment_candidates: dict[str, list[str]]
    linked_candidates: list[str]
    excluded_candidates: list[str]
    valid_combinations: list[dict[str, str]]
    issues: list[dict[str, str]] = field(default_factory=list)
    independent_evidence_units: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _fragment_id(fragment: Any) -> str:
    return str(getattr(fragment, "evidence_id", None) or fragment["fragment_id"])


def _source_id(fragment: Any) -> str:
    return str(getattr(fragment, "source_id", None) or fragment.get("source_id", ""))


def fragment_placements(fragment: Any) -> list[CandidatePlacement]:
    placements: list[CandidatePlacement] = []
    alt = getattr(fragment, "alt", None)
    if alt is None and isinstance(fragment, dict):
        alt = fragment.get("alt", {})
    for key, observation in (alt or {}).items():
        variant = getattr(observation, "variant", None)
        if variant is None and isinstance(observation, dict):
            variant = observation.get("variant")
        if variant is None:
            chrom, pos, ref, _alt = key
            pos = int(pos)
            end = pos + len(str(ref)) - 1
            variant_id = f"{chrom}:{pos}"
        else:
            chrom = getattr(variant, "chrom", None)
            if chrom is None:
                chrom = variant["chrom"]
            pos = int(getattr(variant, "pos", None) or variant["pos"])
            ref = getattr(variant, "ref", None) or variant["ref"]
            end = pos + len(str(ref)) - 1
            variant_id = str(getattr(variant, "vid", None) or variant["variant_id"])
        strand = "F"
        if hasattr(observation, "strand"):
            strand = str(observation.strand)
        elif isinstance(observation, dict):
            strand = str(observation.get("strand", "F"))
        placements.append(CandidatePlacement(variant_id, str(chrom), pos, end, strand))
    return placements


def _index_fragments(fragments: Iterable[Any]) -> tuple[dict[str, Any], list[str]]:
    index: dict[str, Any] = {}
    duplicates: list[str] = []
    for fragment in fragments:
        fid = _fragment_id(fragment)
        if fid in index:
            duplicates.append(fid)
        else:
            index[fid] = fragment
    return index, duplicates


def _orientation_ok(first: CandidatePlacement, second: CandidatePlacement,
                    expected: str) -> bool:
    if expected in ("any", "*"):
        return True
    same = first.strand == second.strand
    if expected == "same":
        return same
    if expected == "opposite":
        return not same
    return False


def _distance_ok(first: CandidatePlacement, second: CandidatePlacement,
                 relation: AssociationRelation) -> bool:
    if first.chrom != second.chrom:
        return False
    distance = second.pos - first.pos
    if relation.min_distance is not None and distance < relation.min_distance:
        return False
    if relation.max_distance is not None and distance > relation.max_distance:
        return False
    return True


def reconcile_associations(
    fragments: Iterable[Any],
    relations: Iterable[AssociationRelation],
) -> list[EvidenceUnitResult]:
    fragment_index, imported_duplicates = _index_fragments(fragments)
    results: list[EvidenceUnitResult] = []
    for relation in relations:
        endpoint_ids = [relation.a, relation.b]
        present = [fid for fid in dict.fromkeys(endpoint_ids) if fid in fragment_index]
        missing = [fid for fid in dict.fromkeys(endpoint_ids) if fid not in fragment_index]
        duplicate_fragment_ids = sorted(set(imported_duplicates) & set(endpoint_ids))

        single = {
            fid: sorted({p.variant_id for p in fragment_placements(fragment_index[fid])})
            for fid in present
        }
        issues: list[dict[str, str]] = []

        if relation.status == "withdrawn":
            results.append(EvidenceUnitResult(
                relation.association_id, relation.molecule_id, relation.relation_id,
                relation.relation_version, "withdrawn", endpoint_ids, present,
                missing, duplicate_fragment_ids, single,
                sorted(set().union(*single.values()) if single else []), [], [],
                [{
                    "code": "association_withdrawn", "fragment_id": "",
                    "message": "relation withdrawn; original ambiguity is restored",
                }], 1 if present else 0))
            continue
        if relation.status == "revised":
            results.append(EvidenceUnitResult(
                relation.association_id, relation.molecule_id, relation.relation_id,
                relation.relation_version, REVISED, endpoint_ids, present,
                missing, duplicate_fragment_ids, single,
                sorted(set().union(*single.values()) if single else []), [], [],
                [{
                    "code": "association_revised", "fragment_id": "",
                    "message": "old relation revised; original ambiguity is retained",
                }], 1 if present else 0))
            continue

        expected_sources = {relation.a: relation.source_id_a, relation.b: relation.source_id_b}
        for fid in present:
            fragment = fragment_index[fid]
            if _source_id(fragment) != expected_sources[fid]:
                issues.append({
                    "code": RELATION_INCONSISTENT,
                    "fragment_id": fid,
                    "message": (
                        f"source {_source_id(fragment)!r} does not match association "
                        f"endpoint source {expected_sources[fid]!r}"),
                })

        if missing:
            status = PENDING
            issues.append({
                "code": PENDING,
                "fragment_id": ",".join(missing),
                "message": "endpoint has not arrived yet; retain unpaired state",
            })
        elif any(not single.get(fid) for fid in present):
            status = MATCH_INSUFFICIENT
            for fid in present:
                if not single.get(fid):
                    issues.append({
                        "code": MATCH_INSUFFICIENT,
                        "fragment_id": fid,
                        "message": "fragment has no usable candidate placement",
                    })
        elif issues:
            status = RELATION_INCONSISTENT
        else:
            placements_a = fragment_placements(fragment_index[relation.a])
            placements_b = fragment_placements(fragment_index[relation.b])
            valid = []
            for first in placements_a:
                for second in placements_b:
                    if (_orientation_ok(first, second, relation.orientation)
                            and _distance_ok(first, second, relation)):
                        valid.append({
                            relation.a: first.variant_id,
                            relation.b: second.variant_id,
                            "chrom": first.chrom,
                            "distance": str(second.pos - first.pos),
                        })
            if not valid:
                status = RELATION_INCONSISTENT
                issues.append({
                    "code": RELATION_INCONSISTENT,
                    "fragment_id": "",
                    "message": "no joint placement satisfies the asserted relation; both endpoints remain visible",
                })
            else:
                status = RESOLVED

        union = sorted(set().union(*single.values()) if single else [])
        if status == RESOLVED:
            linked = sorted({
                candidate
                for combo in valid
                for candidate in (combo[relation.a], combo[relation.b])
            })
            valid_combinations = valid
        else:
            linked = union
            valid_combinations = []

        if duplicate_fragment_ids:
            issues.append({
                "code": "duplicate_import",
                "fragment_id": ",".join(duplicate_fragment_ids),
                "message": "repeated import is collapsed into the same evidence unit",
            })

        results.append(EvidenceUnitResult(
            relation.association_id, relation.molecule_id, relation.relation_id,
            relation.relation_version, status, endpoint_ids, present, missing,
            duplicate_fragment_ids, single, linked,
            sorted(set(union) - set(linked)), valid_combinations, issues,
            1 if present else 0))
    return results


def load_association_relations(path: str | Path | None) -> list[AssociationRelation]:
    if path is None:
        return []
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    records = data.get("relations", data) if isinstance(data, dict) else data
    allowed = AssociationRelation.__dataclass_fields__
    return [
        AssociationRelation(**{k: v for k, v in record.items() if k in allowed})
        for record in records
    ]
