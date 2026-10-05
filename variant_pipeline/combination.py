"""Compatible combination search over multi-source fragment relations.

A researcher declares constraints between association relations: groups
that must appear together (co_required) and groups that must stay
mutually exclusive (mutually_exclusive).  The search enumerates
combinations of relations using (source, fragment) identity only: a
missing endpoint keeps its relation in an explicit unknown state and
is never filled in by a same-named fragment from another source, and a
repeated import of the same (source, fragment) pair never adds support
weight.

The search itself never confirms anything: it produces candidate plans
under two comparable objectives (cover the most relations / use the
fewest independent sources).  A plan only takes effect through an
explicit adoption, which records the selected sources alongside the
plan.  All ordering is canonicalised, so permuting the input relations
or constraints yields an equivalent result.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Iterable

from .store import sha256_text


CO_REQUIRED = "co_required"
MUTUALLY_EXCLUSIVE = "mutually_exclusive"
CONSTRAINT_KINDS = (CO_REQUIRED, MUTUALLY_EXCLUSIVE)

MAX_COVERAGE = "max_coverage"
MIN_SOURCES = "min_sources"
OBJECTIVES = (MAX_COVERAGE, MIN_SOURCES)

# Relations in these states may appear in a plan; anything else (pending
# endpoint, insufficient match, withdrawn, revised) stays visible but is
# not confirmable by the search.
_ELIGIBLE_STATUSES = ("resolved", "molecule_conflict")
_UNKNOWN_STATUSES = ("pending_endpoint",)


@dataclass(frozen=True)
class CombinationConstraint:
    constraint_id: str
    kind: str
    association_ids: tuple

    def __post_init__(self) -> None:
        if self.kind not in CONSTRAINT_KINDS:
            raise ValueError(
                f"unknown constraint kind {self.kind!r} in "
                f"{self.constraint_id!r}; expected one of {CONSTRAINT_KINDS}")
        if len(set(self.association_ids)) != len(self.association_ids):
            raise ValueError(
                f"constraint {self.constraint_id!r} repeats an association id")
        if len(self.association_ids) < 2:
            raise ValueError(
                f"constraint {self.constraint_id!r} needs at least two "
                "association ids")


@dataclass
class RelationOption:
    """One association relation as seen by the combination search."""

    association_id: str
    status: str
    fragment_ids: list
    present_fragment_ids: list
    missing_fragment_ids: list
    sources: list
    relation_id: str | None = None
    relation_version: str = "1"


@dataclass
class CombinationPlan:
    plan_id: str
    objective: str
    association_ids: list
    entity_fragment_ids: list
    entity_count: int
    independent_sources: list
    independent_source_count: int
    covered_relation_count: int
    status: str = "candidate"  # never auto-confirmed by the search

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class CombinationSearchResult:
    status: str  # "ok" | "no_solution"
    plans: dict = field(default_factory=dict)
    comparison: dict = field(default_factory=dict)
    unknown_relations: list = field(default_factory=list)
    blocked_constraints: list = field(default_factory=list)
    conflicts: list = field(default_factory=list)
    constraints: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def relation_options_from_units(units):
    """Build relation options from stored evidence_units records.

    Sources are taken from the canonical source::name endpoint ids, so
    a source re-imported under the same id contributes exactly one entry.
    """
    options = []
    for unit in units:
        present = list(unit.get("present_fragment_ids", []))
        sources = sorted({
            str(fid).split("::", 1)[0] for fid in present if "::" in str(fid)
        })
        options.append(RelationOption(
            association_id=str(unit["association_id"]),
            status=str(unit.get("status", "")),
            fragment_ids=list(unit.get("fragment_ids", [])),
            present_fragment_ids=present,
            missing_fragment_ids=list(unit.get("missing_fragment_ids", [])),
            sources=sources,
            relation_id=unit.get("relation_id"),
            relation_version=str(unit.get("relation_version", "1")),
        ))
    return options


def _canonical_constraints(constraints):
    """Canonical constraint order: input permutation cannot change results."""
    return sorted(
        constraints,
        key=lambda c: (c.kind, tuple(sorted(c.association_ids)),
                       c.constraint_id),
    )


def _satisfies(subset, constraints):
    for constraint in constraints:
        members = set(constraint.association_ids)
        chosen = members & subset
        if constraint.kind == CO_REQUIRED and chosen and chosen != members:
            return False
        if constraint.kind == MUTUALLY_EXCLUSIVE and len(chosen) > 1:
            return False
    return True


def _build_plan(association_ids, options, objective):
    entities = sorted({
        fragment
        for aid in association_ids
        for fragment in options[aid].present_fragment_ids
    })
    sources = sorted({
        source
        for aid in association_ids
        for source in options[aid].sources
    })
    material = json.dumps(list(association_ids), ensure_ascii=False)
    return CombinationPlan(
        plan_id="plan-" + sha256_text(material)[:12],
        objective=objective,
        association_ids=list(association_ids),
        entity_fragment_ids=entities,
        entity_count=len(entities),
        independent_sources=sources,
        independent_source_count=len(sources),
        covered_relation_count=len(association_ids),
    )


def _plan_key(plan_ids, options, objective):
    entities = {
        fragment
        for aid in plan_ids
        for fragment in options[aid].present_fragment_ids
    }
    sources = {
        source
        for aid in plan_ids
        for source in options[aid].sources
    }
    covered = len(plan_ids)
    if objective == MAX_COVERAGE:
        # Most relations covered, then fewest independent sources, then a
        # deterministic lexicographic tie-break on the sorted id tuple.
        return (covered, -len(sources), -len(entities),
                tuple(reversed(plan_ids)))
    # Fewest independent sources, then most relations covered.
    return (-len(sources), covered, -len(entities),
            tuple(reversed(plan_ids)))


def _detect_conflicts(constraints, options):
    """Constraint pairs that cannot both be honoured for their members.

    A co-required group and a mutually-exclusive group overlapping in two
    or more relations demand that those relations both appear together and
    never appear together: no plan can satisfy both, so the conflicting
    constraints are reported with the endpoints of the shared relations.
    """
    conflicts = []
    required = [c for c in constraints if c.kind == CO_REQUIRED]
    exclusive = [c for c in constraints if c.kind == MUTUALLY_EXCLUSIVE]
    for req in required:
        for exc in exclusive:
            shared = sorted(set(req.association_ids) & set(exc.association_ids))
            if len(shared) < 2:
                continue
            endpoints = {
                aid: {
                    "present_fragment_ids": options[aid].present_fragment_ids,
                    "missing_fragment_ids": options[aid].missing_fragment_ids,
                }
                for aid in shared
                if aid in options
            }
            conflicts.append({
                "code": "conflicting_constraints",
                "co_required_constraint": req.constraint_id,
                "mutually_exclusive_constraint": exc.constraint_id,
                "shared_association_ids": shared,
                "endpoints": endpoints,
                "message": (
                    "relations must simultaneously appear together and "
                    "remain mutually exclusive; no combination can satisfy "
                    "both constraints"),
            })
    return sorted(conflicts, key=lambda c: (
        c["co_required_constraint"], c["mutually_exclusive_constraint"]))


def search_combinations(units, constraints):
    """Search compatible relation combinations under both objectives.

    Returns one candidate plan per objective plus a comparison.  The
    search is a pure function of the stored evidence units and the
    declared constraints; it confirms nothing by itself.
    """
    options_list = relation_options_from_units(units)
    options = {opt.association_id: opt for opt in options_list}
    canonical = _canonical_constraints(constraints)

    known_ids = set(options)
    for constraint in canonical:
        unknown = sorted(set(constraint.association_ids) - known_ids)
        if unknown:
            raise ValueError(
                f"constraint {constraint.constraint_id!r} references unknown "
                f"association ids: {', '.join(unknown)}")

    eligible = sorted(
        aid for aid, opt in options.items() if opt.status in _ELIGIBLE_STATUSES)
    unknown_relations = [
        {
            "association_id": opt.association_id,
            "status": opt.status,
            "present_fragment_ids": opt.present_fragment_ids,
            "missing_fragment_ids": opt.missing_fragment_ids,
            "message": (
                "endpoint missing; relation kept as unknown and never "
                "filled by a same-named fragment from another source"),
        }
        for opt in sorted(options_list, key=lambda o: o.association_id)
        if opt.status in _UNKNOWN_STATUSES
    ]
    blocked = []
    for constraint in canonical:
        if constraint.kind != CO_REQUIRED:
            continue
        pending = sorted(
            aid for aid in constraint.association_ids
            if options[aid].status in _UNKNOWN_STATUSES)
        if pending:
            blocked.append({
                "constraint_id": constraint.constraint_id,
                "blocked_by_association_ids": pending,
                "message": (
                    "co-required group cannot be confirmed while a member "
                    "relation has a missing endpoint"),
            })

    conflicts = _detect_conflicts(canonical, options)
    constraint_dicts = [
        {
            "constraint_id": c.constraint_id,
            "kind": c.kind,
            "association_ids": sorted(c.association_ids),
        }
        for c in canonical
    ]
    if conflicts:
        return CombinationSearchResult(
            status="no_solution",
            unknown_relations=unknown_relations,
            blocked_constraints=blocked,
            conflicts=conflicts,
            constraints=constraint_dicts,
        )

    # Exhaustive subset enumeration over the canonically ordered eligible
    # relations; the best plan per objective is unique because ties break
    # on the sorted association id tuple.
    best = {}
    best_key = {}
    count = len(eligible)
    for mask in range(1 << count):
        subset = frozenset(
            eligible[i] for i in range(count) if mask & (1 << i))
        if not _satisfies(subset, canonical):
            continue
        ordered = tuple(sorted(subset))
        for objective in OBJECTIVES:
            # min_sources compares plans that actually cover relations;
            # the empty plan would trivially "use" zero sources.
            if objective == MIN_SOURCES and not ordered:
                continue
            key = _plan_key(ordered, options, objective)
            if objective not in best_key or key > best_key[objective]:
                best_key[objective] = key
                best[objective] = ordered
    # No relation can be covered at all: both objectives fall back to the
    # empty plan so the result stays explicit rather than disappearing.
    best.setdefault(MIN_SOURCES, ())
    best.setdefault(MAX_COVERAGE, ())

    plans = {
        objective: _build_plan(best[objective], options, objective).to_dict()
        for objective in OBJECTIVES
    }
    same = best[MAX_COVERAGE] == best[MIN_SOURCES]
    comparison = {
        "same_plan": same,
        "max_coverage": {
            "plan_id": plans[MAX_COVERAGE]["plan_id"],
            "covered_relation_count": plans[MAX_COVERAGE][
                "covered_relation_count"],
            "independent_source_count": plans[MAX_COVERAGE][
                "independent_source_count"],
        },
        "min_sources": {
            "plan_id": plans[MIN_SOURCES]["plan_id"],
            "covered_relation_count": plans[MIN_SOURCES][
                "covered_relation_count"],
            "independent_source_count": plans[MIN_SOURCES][
                "independent_source_count"],
        },
        "message": (
            "both objectives select the same combination"
            if same else
            "max_coverage covers the most relations; min_sources uses the "
            "fewest independent sources; neither is confirmed until "
            "explicitly adopted"),
    }
    return CombinationSearchResult(
        status="ok",
        plans=plans,
        comparison=comparison,
        unknown_relations=unknown_relations,
        blocked_constraints=blocked,
        constraints=constraint_dicts,
    )


def load_combination_constraints(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    records = data.get("constraints", data) if isinstance(data, dict) else data
    constraints = []
    for record in records:
        constraints.append(CombinationConstraint(
            constraint_id=str(record["constraint_id"]),
            kind=str(record["kind"]),
            association_ids=tuple(str(a) for a in record["association_ids"]),
        ))
    return constraints


def constraints_fingerprint(constraints):
    canonical = _canonical_constraints(constraints)
    material = json.dumps(
        [
            {
                "constraint_id": c.constraint_id,
                "kind": c.kind,
                "association_ids": sorted(c.association_ids),
            }
            for c in canonical
        ],
        sort_keys=True, ensure_ascii=False)
    return sha256_text(material)

COMBINATIONS_DIRNAME = "combinations"


def combination_path(store, version: str) -> Path:
    safe = version.replace("/", "_")
    return store.root / COMBINATIONS_DIRNAME / f"{safe}.json"


def load_combination_adoption(store, version: str) -> dict | None:
    path = combination_path(store, version)
    if not path.is_file():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("status") == "adopted" else None


def adopt_combination(
    store,
    version: str,
    constraints: list,
    plan_id: str,
    rationale: str | None = None,
) -> dict:
    """Explicitly adopt one candidate plan; the base version is untouched.

    The search is recomputed from the stored evidence units so the adopted
    plan is always the deterministic result of the recorded constraints;
    the record keeps the selected sources alongside the plan.  A previous
    adoption is archived into history, never overwritten.
    """
    info = store.verify(version)
    payload = json.loads(
        (info.path / "candidates.json").read_text(encoding="utf-8"))
    result = search_combinations(payload.get("evidence_units", []), constraints)
    if result.status != "ok":
        raise ValueError(
            f"no adoptable combination for {version}: conflicting "
            "constraints; see search conflicts for endpoints")
    plan = next(
        (p for p in result.plans.values() if p["plan_id"] == plan_id), None)
    if plan is None:
        known = sorted(p["plan_id"] for p in result.plans.values())
        raise ValueError(
            f"unknown plan {plan_id!r} for {version}; candidate plans: "
            f"{', '.join(known)}")

    target = combination_path(store, version)
    existing = (json.loads(target.read_text(encoding="utf-8"))
                if target.is_file() else None)
    history = list(existing.get("history", [])) if existing else []
    if existing and existing.get("status") == "adopted":
        history.append({k: existing.get(k) for k in (
            "adoption_id", "plan_id", "objective", "association_ids",
            "selected_sources", "constraints_fingerprint", "rationale",
            "created_at")})

    import os
    import tempfile
    from datetime import datetime, timezone
    record = {
        "adoption_id": "adopt-" + sha256_text(
            f"{version}:{plan_id}:{len(history)}")[:12],
        "base_version": version,
        "status": "adopted",
        "does_not_modify_base_version": True,
        "plan_id": plan["plan_id"],
        "objective": plan["objective"],
        "association_ids": plan["association_ids"],
        "entity_fragment_ids": plan["entity_fragment_ids"],
        "entity_count": plan["entity_count"],
        "selected_sources": plan["independent_sources"],
        "independent_source_count": plan["independent_source_count"],
        "covered_relation_count": plan["covered_relation_count"],
        "constraints_fingerprint": constraints_fingerprint(constraints),
        "rationale": rationale,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "supersedes": existing.get("adoption_id") if existing else None,
        "history": history,
    }
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".combination-", dir=str(target.parent))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, target)
    return record
