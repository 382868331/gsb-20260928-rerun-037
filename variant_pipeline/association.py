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
MOLECULE_CONFLICT = "molecule_conflict"

INTRA_MOLECULE_CONFLICT = "intra_molecule_conflict"


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
    conflicting_candidates: list[str] = field(default_factory=list)
    conflicting_combinations: list[dict[str, Any]] = field(default_factory=list)
    excluded_combinations: list[dict[str, str]] = field(default_factory=list)
    issues: list[dict[str, str]] = field(default_factory=list)
    independent_evidence_units: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _source_id(fragment: Any) -> str:
    return str(getattr(fragment, "source_id", None) or fragment.get("source_id", ""))


def _local_id(fragment: Any) -> str:
    """Fragment identifier within its source (display name, no source prefix)."""
    qname = getattr(fragment, "qname", None)
    if qname is not None:
        return str(qname)
    raw = fragment.get("fragment_id", "") if isinstance(fragment, dict) else ""
    raw = str(raw)
    source = _source_id(fragment)
    prefix = f"{source}::"
    if source and raw.startswith(prefix):
        return raw[len(prefix):]
    return raw


def _qualified_id(fragment: Any) -> str:
    """Source-qualified identity: same display name from different sources
    stays distinct; only the same (source, id) pair collapses."""
    evidence_id = getattr(fragment, "evidence_id", None)
    if evidence_id:
        return str(evidence_id)
    source = _source_id(fragment)
    local = _local_id(fragment)
    return f"{source}::{local}" if source else local


def _qualify_endpoint(endpoint: str, source: str) -> str:
    """Source-qualified form of a relation endpoint for diagnostics."""
    if not source or "::" in endpoint:
        return endpoint
    return f"{source}::{endpoint}"


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


def _index_fragments(
    fragments: Iterable[Any],
) -> tuple[dict[str, Any], dict[tuple[str, str], str], list[str]]:
    """Index fragments by qualified identity and by (source, local id).

    The same display name imported from different sources is kept as
    distinct fragments; only a repeated (source, id) pair is a duplicate
    import, and the first occurrence wins deterministically.
    """
    by_qualified: dict[str, Any] = {}
    by_pair: dict[tuple[str, str], str] = {}
    duplicates: list[str] = []
    for fragment in fragments:
        qualified = _qualified_id(fragment)
        if qualified in by_qualified:
            duplicates.append(qualified)
            continue
        by_qualified[qualified] = fragment
        by_pair[(_source_id(fragment), _local_id(fragment))] = qualified
    return by_qualified, by_pair, duplicates


def _resolve_endpoint(
    endpoint: str,
    source: str,
    by_qualified: dict[str, Any],
    by_pair: dict[tuple[str, str], str],
) -> str | None:
    """Locate an endpoint by source and fragment identifier jointly.

    Accepts either a source-qualified evidence id or a bare display name
    combined with the endpoint's declared source; never falls back to a
    display-name-only match across sources.
    """
    if endpoint in by_qualified:
        return endpoint
    return by_pair.get((source, endpoint))


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


def _placements_overlap(first: CandidatePlacement,
                        second: CandidatePlacement) -> bool:
    """Whether two allele events claim the same reference interval."""
    return (
        first.chrom == second.chrom
        and first.pos <= second.end
        and second.pos <= first.end
    )


def reconcile_associations(
    fragments: Iterable[Any],
    relations: Iterable[AssociationRelation],
) -> list[EvidenceUnitResult]:
    by_qualified, by_pair, imported_duplicates = _index_fragments(fragments)
    results: list[EvidenceUnitResult] = []
    for relation in relations:
        endpoint_specs = [
            (relation.a, relation.source_id_a),
            (relation.b, relation.source_id_b),
        ]
        resolved: dict[str, str | None] = {}
        for endpoint, source in endpoint_specs:
            key = f"{source}::{endpoint}"
            if key not in resolved:
                resolved[key] = _resolve_endpoint(
                    endpoint, source, by_qualified, by_pair)
        # Canonical endpoint identity is source-qualified, so same-named
        # fragments from different sources stay distinguishable everywhere
        # (status, combinations, diagnostics).
        endpoint_ids = [
            resolved[f"{source}::{endpoint}"] or _qualify_endpoint(endpoint, source)
            for endpoint, source in endpoint_specs
        ]
        present = [
            fid for fid in dict.fromkeys(endpoint_ids) if fid in by_qualified]
        missing = [
            fid for fid in dict.fromkeys(endpoint_ids) if fid not in by_qualified]
        duplicate_fragment_ids = sorted(set(imported_duplicates) & set(present))

        single = {
            fid: sorted({p.variant_id for p in fragment_placements(by_qualified[fid])})
            for fid in present
        }
        issues: list[dict[str, str]] = []
        valid: list[dict[str, str]] = []
        relation_excluded: list[dict[str, str]] = []
        conflicting_combinations: list[dict[str, Any]] = []
        conflicting_variant_ids: set[str] = set()

        if relation.status == "withdrawn":
            results.append(EvidenceUnitResult(
                relation.association_id, relation.molecule_id, relation.relation_id,
                relation.relation_version, "withdrawn", endpoint_ids, present,
                missing, duplicate_fragment_ids, single,
                sorted(set().union(*single.values()) if single else []), [], [],
                issues=[{
                    "code": "association_withdrawn", "fragment_id": "",
                    "message": "relation withdrawn; original ambiguity is restored",
                }],
                independent_evidence_units=1 if present else 0))
            continue
        if relation.status == "revised":
            results.append(EvidenceUnitResult(
                relation.association_id, relation.molecule_id, relation.relation_id,
                relation.relation_version, REVISED, endpoint_ids, present,
                missing, duplicate_fragment_ids, single,
                sorted(set().union(*single.values()) if single else []), [], [],
                issues=[{
                    "code": "association_revised", "fragment_id": "",
                    "message": "old relation revised; original ambiguity is retained",
                }],
                independent_evidence_units=1 if present else 0))
            continue

        for (endpoint, expected_source), qualified in zip(
                endpoint_specs, endpoint_ids):
            if qualified not in by_qualified:
                continue
            fragment = by_qualified[qualified]
            if _source_id(fragment) != expected_source:
                issues.append({
                    "code": RELATION_INCONSISTENT,
                    "fragment_id": qualified,
                    "message": (
                        f"source {_source_id(fragment)!r} does not match association "
                        f"endpoint source {expected_source!r}"),
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
            placements_a = fragment_placements(by_qualified[endpoint_ids[0]])
            placements_b = fragment_placements(by_qualified[endpoint_ids[1]])
            all_pairs: list[tuple[CandidatePlacement, CandidatePlacement,
                                  dict[str, str]]] = []
            for first in placements_a:
                for second in placements_b:
                    pair = {
                        endpoint_ids[0]: first.variant_id,
                        endpoint_ids[1]: second.variant_id,
                        "chrom": first.chrom,
                        "distance": str(second.pos - first.pos),
                    }
                    # Physical allele disagreement inside one molecule cannot
                    # be resolved by the relation: both alleles stay on record
                    # and the conflicting pair is retained as audit evidence.
                    if (_placements_overlap(first, second)
                            and first.variant_id != second.variant_id):
                        conflicting_variant_ids.add(first.variant_id)
                        conflicting_variant_ids.add(second.variant_id)
                        conflict_entry = {
                            "code": INTRA_MOLECULE_CONFLICT,
                            "a_variant_id": first.variant_id,
                            "b_variant_id": second.variant_id,
                            "chrom": first.chrom,
                            "pos": max(first.pos, second.pos),
                            "end": min(first.end, second.end),
                            "fragment_id": "",
                            "message": (
                                "endpoints assert different alleles over the "
                                "same reference interval in one molecule; "
                                "neither allele is dropped or auto-selected"),
                        }
                        conflicting_combinations.append(conflict_entry)
                    else:
                        all_pairs.append((first, second, pair))
            # An allele contradicted at its own site cannot be carried into a
            # joint interpretation via the other endpoint's placements: only
            # combinations free of conflicting alleles are retested.
            for first, second, pair in all_pairs:
                if (first.variant_id in conflicting_variant_ids
                        or second.variant_id in conflicting_variant_ids):
                    continue
                if (_orientation_ok(first, second, relation.orientation)
                        and _distance_ok(first, second, relation)):
                    valid.append(pair)
                else:
                    relation_excluded.append(pair)
            if not valid and not conflicting_variant_ids:
                status = RELATION_INCONSISTENT
                issues.append({
                    "code": RELATION_INCONSISTENT,
                    "fragment_id": "",
                    "message": "no joint placement satisfies the asserted relation; both endpoints remain visible",
                })
            else:
                status = RESOLVED
            if conflicting_combinations:
                status = MOLECULE_CONFLICT
                issues.append({
                    "code": INTRA_MOLECULE_CONFLICT,
                    "fragment_id": "",
                    "message": (
                        "intra-molecule allele conflict at an overlapping "
                        "site; both candidate alleles are retained and "
                        "abstain at that site; the molecule still counts once"),
                })

        union = sorted(set().union(*single.values()) if single else [])
        if status in (RESOLVED, MOLECULE_CONFLICT):
            linked = sorted({
                candidate
                for combo in valid
                for candidate in (combo[endpoint_ids[0]], combo[endpoint_ids[1]])
            })
            valid_combinations = valid
            relation_excluded_sorted = sorted(
                relation_excluded,
                key=lambda combo: (
                    str(combo.get(endpoint_ids[0], "")),
                    str(combo.get(endpoint_ids[1], "")),
                    str(combo.get("chrom", "")),
                    str(combo.get("distance", ""))))
        else:
            linked = union
            valid_combinations = []
            relation_excluded_sorted = []

        conflicting_candidates = sorted(conflicting_variant_ids) if single else []

        if duplicate_fragment_ids:
            issues.append({
                "code": "duplicate_import",
                "fragment_id": ",".join(duplicate_fragment_ids),
                "message": "repeated import is collapsed into the same evidence unit",
            })

        results.append(EvidenceUnitResult(
            association_id=relation.association_id,
            molecule_id=relation.molecule_id,
            relation_id=relation.relation_id,
            relation_version=relation.relation_version,
            status=status,
            fragment_ids=endpoint_ids,
            present_fragment_ids=present,
            missing_fragment_ids=missing,
            duplicate_fragment_ids=duplicate_fragment_ids,
            single_fragment_candidates=single,
            linked_candidates=linked,
            excluded_candidates=sorted(
                set(union) - set(linked) - set(conflicting_candidates)),
            valid_combinations=valid_combinations,
            conflicting_candidates=conflicting_candidates,
            conflicting_combinations=(
                conflicting_combinations if single else []),
            excluded_combinations=relation_excluded_sorted,
            issues=issues,
            independent_evidence_units=1 if present else 0))
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
