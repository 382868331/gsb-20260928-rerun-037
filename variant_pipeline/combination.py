"""Compatible combination search across multiple sourced relations.

Researchers can require that several explicit fragment relations hold
together (``required``) or be mutually exclusive (``forbidden``).  This
layer searches for joint assignments of one valid endpoint->variant
combination per relation that honour those constraints.

Design rules:

* endpoints are matched by (source, fragment identity) together; a same
  named fragment from another source can never substitute for a missing one;
* a fragment shared by several relations is one entity and must assert the
  same candidate everywhere, otherwise the assignments are incompatible;
* repeated imports of one (source, fragment) pair were already collapsed by
  the association layer, so they never add a source or extra support;
* an endpoint that has not arrived keeps the relation *unknown* instead of
  silently dropping it or substituting another fragment;
* searching never confirms a relation; adoption is a separate, explicit act.

Results do not depend on relation, fragment or constraint input ordering:
every collection is canonicalized before enumeration and reporting.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Iterable

from .association import (
    PENDING,
    AssociationRelation,
    reconcile_associations,
)

FEASIBLE = "feasible"
INFEASIBLE = "infeasible"
UNKNOWN = "unknown"  # feasible so far, but a required relation misses an end

REQUIRED = "required"      # all listed relations must co-occur
FORBIDDEN = "forbidden"    # listed relations are mutually exclusive

SEARCHABLE_STATUSES = {"resolved", "molecule_conflict"}


@dataclass(frozen=True)
class CombinationConstraint:
    constraint_id: str
    kind: str  # required | forbidden
    relations: tuple[str, ...]

    @staticmethod
    def from_record(record: dict[str, Any]) -> "CombinationConstraint":
        kind = str(record.get("kind", "")).strip().lower()
        if kind not in (REQUIRED, FORBIDDEN):
            raise ValueError(
                f"combination constraint {record.get('constraint_id')!r} "
                f"has unsupported kind {kind!r}; expected required/forbidden")
        relations = record.get("relations") or record.get("relation_ids")
        if not isinstance(relations, list) or not relations:
            raise ValueError(
                f"combination constraint {record.get('constraint_id')!r} "
                "must list at least one relation id")
        return CombinationConstraint(
            constraint_id=str(record["constraint_id"]),
            kind=kind,
            relations=tuple(str(r) for r in relations),
        )


@dataclass
class CombinationSolution:
    assignments: dict[str, dict[str, str]]  # association_id -> endpoint -> variant
    covered_relations: list[str]
    fragment_entities: list[str]
    independent_sources: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ComboSearchResult:
    status: str
    solutions: list[CombinationSolution]
    unknown_relations: list[str] = field(default_factory=list)
    unsatisfiable_relations: list[str] = field(default_factory=list)
    constraint_reports: list[dict[str, Any]] = field(default_factory=list)
    conflicts: list[dict[str, Any]] = field(default_factory=list)
    adopted: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["solutions"] = [s.to_dict() for s in self.solutions]
        return data


def load_combination_constraints(
    path: str | Path | None,
) -> list[CombinationConstraint]:
    if path is None:
        return []
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    records = data if isinstance(data, list) else (
        data.get("combination_constraints") or data.get("constraints"))
    if not records:
        return []
    return [CombinationConstraint.from_record(r) for r in records]


def constraints_from_association_doc(
    path: str | Path | None,
) -> list[CombinationConstraint]:
    """Read optional inline constraints from the relations document."""
    if path is None:
        return []
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(data, dict):
        records = data.get("combination_constraints") or data.get("constraints")
        if records:
            return [CombinationConstraint.from_record(r) for r in records]
    return []


def _source_of(endpoint_id: str) -> str:
    return endpoint_id.partition("::")[0] if "::" in endpoint_id else ""


def _searchable_units(units: list[Any]) -> dict[str, Any]:
    return {
        unit.association_id: unit
        for unit in units
        if unit.status in SEARCHABLE_STATUSES and unit.valid_combinations
    }


def search_combinations(
    fragments: Iterable[Any],
    relations: Iterable[AssociationRelation],
    constraints: Iterable[CombinationConstraint] = (),
) -> ComboSearchResult:
    units = reconcile_associations(fragments, list(relations))
    units_by_id = {unit.association_id: unit for unit in units}
    searchable = _searchable_units(units)

    constraints = sorted(constraints, key=lambda c: c.constraint_id)
    required_ids = sorted({
        rid for c in constraints if c.kind == REQUIRED for rid in c.relations
    })
    forbidden_groups = sorted(
        (sorted(set(c.relations)) for c in constraints if c.kind == FORBIDDEN),
        key=lambda group: tuple(group),
    )

    referenced = sorted({rid for c in constraints for rid in c.relations})

    unknown_relations: set[str] = set()
    unsatisfiable_relations: set[str] = set()
    pending_ids = {
        unit.association_id
        for unit in units
        if unit.status == PENDING or unit.missing_fragment_ids
    }
    for rid in sorted(set(referenced) | set(searchable) | pending_ids):
        unit = units_by_id.get(rid)
        if unit is None:
            # A constraint names a relation that was never supplied: there is
            # no way to satisfy it and no same-named stand-in may be used.
            unsatisfiable_relations.add(rid)
        elif unit.status == PENDING or unit.missing_fragment_ids:
            unknown_relations.add(rid)
        elif rid not in searchable and rid in referenced:
            unsatisfiable_relations.add(rid)

    conflicts: list[dict[str, Any]] = []
    hard_required_missing = sorted(set(required_ids) & unsatisfiable_relations)
    if hard_required_missing:
        for rid in hard_required_missing:
            unit = units_by_id.get(rid)
            endpoints = unit.fragment_ids if unit else []
            for c in constraints:
                if c.kind == REQUIRED and rid in c.relations:
                    conflicts.append({
                        "code": "required_relation_unavailable",
                        "constraint_id": c.constraint_id,
                        "relation_id": rid,
                        "endpoints": endpoints,
                        "message": (
                            "required relation has no joint placement and no "
                            "same-named fragment may substitute"),
                    })
        return ComboSearchResult(
            status=INFEASIBLE, solutions=[],
            unknown_relations=sorted(unknown_relations),
            unsatisfiable_relations=sorted(unsatisfiable_relations),
            constraint_reports=[],
            conflicts=conflicts)

    relation_ids = sorted(searchable)
    domains = {rid: searchable[rid].valid_combinations for rid in relation_ids}
    required_set = set(required_ids)

    found: list[dict[str, dict[str, str]]] = []

    def compatible(assignment, rid, combo):
        # Shared endpoint consistency: one fragment asserts one candidate.
        for other_combo in assignment.values():
            for endpoint_id, variant_id in combo.items():
                if endpoint_id in other_combo and                         other_combo[endpoint_id] != variant_id:
                    return False
        # Mutual exclusion: at most one relation of each forbidden group.
        for group in forbidden_groups:
            if rid in group and any(other in group for other in assignment):
                return False
        return True

    def visit(index, assignment):
        if index == len(relation_ids):
            found.append(dict(assignment))
            return
        rid = relation_ids[index]
        meta_keys = {"chrom", "distance"}
        combos = sorted(
            ({k: v for k, v in combo.items() if k not in meta_keys}
             for combo in domains[rid]),
            key=lambda combo: tuple(sorted(combo.items())),
        )
        for combo in combos:
            if compatible(assignment, rid, combo):
                assignment[rid] = combo
                visit(index + 1, assignment)
                del assignment[rid]
        if rid not in required_set:  # non-required relations may be skipped
            visit(index + 1, assignment)

    visit(0, {})

    def summarize(assignment):
        entities = sorted({ep for combo in assignment.values() for ep in combo})
        sources = sorted({_source_of(ep) for ep in entities})
        return CombinationSolution(
            assignments={rid: dict(assignment[rid]) for rid in sorted(assignment)},
            covered_relations=sorted(assignment),
            fragment_entities=entities,
            independent_sources=sources,
        )

    solutions = [summarize(a) for a in found]
    solutions.sort(key=lambda s: (
        -len(s.covered_relations),
        len(s.independent_sources),
        len(s.fragment_entities),
        json.dumps(s.assignments, sort_keys=True, ensure_ascii=False),
    ))

    required_unknown = bool(set(required_ids) & unknown_relations)
    if not solutions:
        conflicts.extend(_forbidden_conflicts(constraints, units_by_id, searchable))
        status = INFEASIBLE
    elif required_unknown:
        status = UNKNOWN
    else:
        status = FEASIBLE

    constraint_reports = _constraint_reports(
        constraints, solutions, units_by_id, searchable,
        unknown_relations, unsatisfiable_relations)

    return ComboSearchResult(
        status=status,
        solutions=solutions,
        unknown_relations=sorted(unknown_relations),
        unsatisfiable_relations=sorted(unsatisfiable_relations),
        constraint_reports=constraint_reports,
        conflicts=conflicts)


def _forbidden_conflicts(constraints, units_by_id, searchable):
    reports = []
    for c in constraints:
        if c.kind != FORBIDDEN:
            continue
        present = [rid for rid in sorted(set(c.relations)) if rid in searchable]
        if len(present) >= 2:
            endpoints = sorted({
                endpoint_id
                for rid in present
                for endpoint_id in units_by_id[rid].fragment_ids
            })
            reports.append({
                "code": "forbidden_relations_cooccur",
                "constraint_id": c.constraint_id,
                "relation_id": ",".join(present),
                "endpoints": endpoints,
                "message": (
                    "mutually exclusive relations can each be placed, so no "
                    "joint combination satisfies the constraint"),
            })
    return reports


def _constraint_reports(constraints, solutions, units_by_id, searchable,
                        unknown_relations, unsatisfiable_relations):
    best = solutions[0] if solutions else None
    covered = set(best.covered_relations) if best else set()
    reports = []
    for c in constraints:
        rids = sorted(set(c.relations))
        endpoints = sorted({
            endpoint
            for rid in rids if rid in units_by_id
            for endpoint in units_by_id[rid].fragment_ids
        })
        if any(rid in unsatisfiable_relations for rid in rids):
            state = "unsatisfiable"
        elif any(rid in unknown_relations for rid in rids):
            state = "unknown"
        elif c.kind == REQUIRED:
            state = "satisfied" if all(rid in covered for rid in rids) else "violated"
        else:
            co_occurring = sum(1 for rid in rids if rid in covered)
            state = "satisfied" if co_occurring <= 1 else "violated"
        reports.append({
            "constraint_id": c.constraint_id,
            "kind": c.kind,
            "relations": rids,
            "endpoints": endpoints,
            "state": state,
        })
    reports.sort(key=lambda r: r["constraint_id"])
    return reports


def adopt_solution(
    result: ComboSearchResult,
    selected_index: int = 0,
    rationale: str | None = None,
) -> dict[str, Any]:
    """Explicitly adopt one searched solution without mutating the run.

    Candidate search never confirms anything on its own; the chosen sources
    and endpoint assignments are recorded only when a researcher adopts a
    ranked solution.
    """
    if result.status == INFEASIBLE:
        raise ValueError("cannot adopt a combination: search is infeasible")
    if not 0 <= selected_index < len(result.solutions):
        raise ValueError(
            f"solution index {selected_index} out of range "
            f"(0..{len(result.solutions) - 1})")
    solution = result.solutions[selected_index]
    return {
        "selected_index": selected_index,
        "status": result.status,
        "rationale": rationale,
        "assignments": solution.assignments,
        "covered_relations": solution.covered_relations,
        "fragment_entities": solution.fragment_entities,
        "independent_sources": solution.independent_sources,
    }
