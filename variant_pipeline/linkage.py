"""Linked fragments from one original molecule as a single evidence unit.

A linkage record (--links links.json) declares that several fragments
came from the same original molecule. Linked fragments form **one**
evidence unit: they never count as multiple independent supports, and the
combination of their individual matches can exclude interpretations that
each fragment alone would leave open.

Linkage information may arrive late or be revised:

- Only one end observed -> the unit stays in 'pending' state; the
  observed members still form a single explainable unit.
- Same qname but a different source / molecule identity -> never paired
  automatically; the link is reported as 'ambiguous_identity'.
- Declared direction / distance / source relations contradicting the
  observed alignments -> 'inconsistent'; both ends are kept visible as
  fragment-level evidence and the contradiction is reported, never
  silently dropped. Reasons distinguish 'fragment_matching_insufficient'
  (a member has no usable alignment) from
  'linkage_evidence_inconsistent' (the declared relation is violated).
- Re-importing the same molecule (duplicate molecule_id) merges into
  one unit and is recorded as 'duplicate_import'.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

ORIENTATIONS = ("FR", "same", "any")

# Link/unit states
STATE_COMPLETE = "complete"
STATE_PENDING = "pending"
STATE_INCONSISTENT = "inconsistent"
STATE_AMBIGUOUS_IDENTITY = "ambiguous_identity"
STATE_SINGLETON = "singleton"

# Reason kinds
REASON_FRAGMENT_INSUFFICIENT = "fragment_matching_insufficient"
REASON_LINKAGE_INCONSISTENT = "linkage_evidence_inconsistent"
REASON_DUPLICATE_IMPORT = "duplicate_import"
REASON_IDENTITY_AMBIGUOUS = "identity_ambiguous"


@dataclass(frozen=True)
class LinkSpec:
    """One declared linkage between fragments of the same molecule."""

    link_id: str
    molecule_id: str
    members: tuple[str, ...]
    source: str | None = None
    expected_orientation: str = "any"
    max_distance: int | None = None
    same_chrom: bool = False
    origin: str | None = None  # links file this record came from


def load_links(path: str | Path) -> list[LinkSpec]:
    """Parse a links JSON file: {"links": [{...}, ...]}."""
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    records = data.get("links") if isinstance(data, dict) else data
    if not isinstance(records, list):
        raise ValueError(f"links file {path} must contain a 'links' list")
    specs: list[LinkSpec] = []
    for index, item in enumerate(records, start=1):
        if not isinstance(item, dict):
            raise ValueError(f"links file {path}: entry {index} is not an object")
        members = item.get("members")
        if not isinstance(members, list) or not members:
            raise ValueError(
                f"links file {path}: entry {index} needs a non-empty 'members' list")
        expected = item.get("expected", {}) or {}
        orientation = expected.get("orientation", "any")
        if orientation not in ORIENTATIONS:
            raise ValueError(
                f"links file {path}: unknown orientation {orientation!r} "
                f"(expected one of {', '.join(ORIENTATIONS)})")
        max_distance = expected.get("max_distance")
        specs.append(LinkSpec(
            link_id=str(item.get("link_id") or f"link-{index:04d}"),
            molecule_id=str(item.get("molecule_id") or item.get("link_id")
                            or f"mol-{index:04d}"),
            members=tuple(str(m) for m in members),
            source=item.get("source"),
            expected_orientation=orientation,
            max_distance=int(max_distance) if max_distance is not None else None,
            same_chrom=bool(expected.get("same_chrom", False)),
            origin=str(path),
        ))
    return specs


@dataclass
class FragmentView:
    """Minimal per-fragment linkage facts (usable from live fragments or
    from the detail records stored in a published payload)."""

    fragment_id: str
    qname: str
    source_id: str
    source_name: str
    strand: str
    chrom: str | None
    spans: tuple[tuple[int, int], ...]
    alt_variant_ids: tuple[str, ...]

    @staticmethod
    def from_detail(detail: dict) -> "FragmentView":
        return FragmentView(
            fragment_id=detail["fragment_id"],
            qname=detail["qname"],
            source_id=detail.get("source_id") or "",
            source_name=detail.get("source_name") or "",
            strand=detail.get("strand") or "F",
            chrom=detail.get("chrom"),
            spans=tuple(tuple(s) for s in detail.get("spans", [])),
            alt_variant_ids=tuple(detail.get("alt_variant_ids", [])),
        )


@dataclass
class LinkUnit:
    unit_id: str
    link_id: str | None
    molecule_id: str | None
    state: str
    members: list[FragmentView] = field(default_factory=list)
    missing_members: list[str] = field(default_factory=list)
    reasons: list[dict] = field(default_factory=list)
    origin: str | None = None

    def to_dict(self) -> dict:
        return {
            "unit_id": self.unit_id,
            "link_id": self.link_id,
            "molecule_id": self.molecule_id,
            "state": self.state,
            "members": [
                {
                    "fragment_id": m.fragment_id,
                    "qname": m.qname,
                    "source_id": m.source_id,
                    "source_name": m.source_name,
                    "strand": m.strand,
                    "chrom": m.chrom,
                    "spans": [list(s) for s in m.spans],
                    "alt_variant_ids": list(m.alt_variant_ids),
                }
                for m in self.members
            ],
            "missing_members": list(self.missing_members),
            "reasons": list(self.reasons),
            "origin": self.origin,
        }


def _source_matches(spec: LinkSpec, view: FragmentView) -> bool:
    if spec.source is None:
        return True
    return spec.source in (view.source_id, view.source_name)


def _span_distance(a: FragmentView, b: FragmentView) -> int | None:
    """Gap in bp between the nearest spans of two same-chrom fragments."""
    if a.chrom is None or b.chrom is None or a.chrom != b.chrom:
        return None
    best: int | None = None
    for start_a, end_a in a.spans:
        for start_b, end_b in b.spans:
            gap = max(0, max(start_a, start_b) - min(end_a, end_b) - 1)
            if best is None or gap < best:
                best = gap
    return best


def _check_constraints(spec: LinkSpec, members: list[FragmentView]) -> list[dict]:
    reasons: list[dict] = []
    for member in members:
        if member.chrom is None or not member.spans:
            reasons.append({
                "kind": REASON_FRAGMENT_INSUFFICIENT,
                "detail": f"member {member.fragment_id} has no usable alignment span",
            })
    if spec.same_chrom and len({m.chrom for m in members}) > 1:
        reasons.append({
            "kind": REASON_LINKAGE_INCONSISTENT,
            "detail": "members align to different references but same_chrom "
                      "was declared",
        })
    if spec.max_distance is not None and len(members) >= 2:
        ordered = sorted(
            members, key=lambda m: m.spans[0][0] if m.spans else 0)
        for first, second in zip(ordered, ordered[1:]):
            gap = _span_distance(first, second)
            if gap is None:
                reasons.append({
                    "kind": REASON_LINKAGE_INCONSISTENT,
                    "detail": "cannot measure distance across references for "
                              f"{first.fragment_id} and {second.fragment_id}",
                })
            elif gap > spec.max_distance:
                reasons.append({
                    "kind": REASON_LINKAGE_INCONSISTENT,
                    "detail": f"observed distance {gap} bp between "
                              f"{first.fragment_id} and {second.fragment_id} "
                              f"exceeds declared max_distance "
                              f"{spec.max_distance} bp",
                })
    if spec.expected_orientation != "any" and len(members) >= 2:
        ordered = sorted(
            members, key=lambda m: m.spans[0][0] if m.spans else 0)
        strands = [m.strand for m in ordered]
        if spec.expected_orientation == "same" and len(set(strands)) > 1:
            reasons.append({
                "kind": REASON_LINKAGE_INCONSISTENT,
                "detail": "declared orientation 'same' but observed strands "
                          f"{'/'.join(strands)}",
            })
        elif (spec.expected_orientation == "FR" and len(strands) == 2
                and strands != ["F", "R"]):
            reasons.append({
                "kind": REASON_LINKAGE_INCONSISTENT,
                "detail": "declared orientation 'FR' but observed strands "
                          f"{'/'.join(strands)} in reference order",
            })
    return reasons


def resolve_links(
    specs: list[LinkSpec],
    views: list[FragmentView],
) -> tuple[list[LinkUnit], list[list[FragmentView]]]:
    """Resolve declared links against observed fragments.

    Returns (all units incl. singletons, counting units): counting units
    are the evidence units used for molecule-level support. Complete and
    pending links each contribute one unit; inconsistent or
    identity-ambiguous links degrade to per-fragment singletons so neither
    end is silently merged or dropped.
    """
    by_qname: dict[str, list[FragmentView]] = {}
    for view in views:
        by_qname.setdefault(view.qname, []).append(view)

    # Merge duplicate declarations of the same molecule into one link.
    merged: dict[str, LinkSpec] = {}
    duplicate_notes: dict[str, list[str]] = {}
    order: list[str] = []
    for spec in specs:
        key = spec.molecule_id
        if key not in merged:
            merged[key] = spec
            duplicate_notes[key] = []
            order.append(key)
        else:
            base = merged[key]
            extra = [m for m in spec.members if m not in base.members]
            merged[key] = LinkSpec(
                link_id=base.link_id,
                molecule_id=base.molecule_id,
                members=base.members + tuple(extra),
                source=base.source or spec.source,
                expected_orientation=base.expected_orientation,
                max_distance=base.max_distance,
                same_chrom=base.same_chrom,
                origin=base.origin,
            )
            duplicate_notes[key].append(spec.link_id)

    assigned: set[str] = set()
    units: list[LinkUnit] = []
    counting: list[list[FragmentView]] = []

    for key in order:
        spec = merged[key]
        unit = LinkUnit(
            unit_id=f"unit-{spec.molecule_id}",
            link_id=spec.link_id,
            molecule_id=spec.molecule_id,
            state=STATE_PENDING,
            origin=spec.origin,
        )
        for dup in duplicate_notes[key]:
            unit.reasons.append({
                "kind": REASON_DUPLICATE_IMPORT,
                "detail": f"molecule {spec.molecule_id} declared again by "
                          f"link {dup}; merged into one evidence unit",
            })
        identity_ambiguous = False
        for qname in spec.members:
            pool = [
                v for v in by_qname.get(qname, [])
                if v.fragment_id not in assigned
            ]
            if (spec.source is None
                    and len({v.source_id for v in pool}) > 1):
                # Same name from multiple sources and no declared
                # source: never auto-pair.
                identity_ambiguous = True
                unit.reasons.append({
                    "kind": REASON_IDENTITY_AMBIGUOUS,
                    "detail": f"member {qname!r} matches fragments from "
                              f"multiple sources "
                              f"{sorted({v.source_name for v in pool})}; "
                              "declare 'source' to pair explicitly",
                })
                unit.missing_members.append(qname)
                continue
            eligible = [v for v in pool if _source_matches(spec, v)]
            if not eligible:
                if pool and spec.source is not None:
                    unit.reasons.append({
                        "kind": REASON_LINKAGE_INCONSISTENT,
                        "detail": f"member {qname!r} observed only from "
                                  f"source(s) "
                                  f"{sorted({v.source_name for v in pool})}, "
                                  f"not declared source {spec.source!r}",
                    })
                    unit.missing_members.append(qname)
                    continue
                unit.reasons.append({
                    "kind": REASON_FRAGMENT_INSUFFICIENT,
                    "detail": f"member {qname!r} has no observed fragment "
                              "(pending or missing)",
                })
                unit.missing_members.append(qname)
                continue
            view = eligible[0]
            assigned.add(view.fragment_id)
            unit.members.append(view)

        if identity_ambiguous:
            unit.state = STATE_AMBIGUOUS_IDENTITY
        elif unit.missing_members:
            unit.state = STATE_PENDING
        else:
            violations = _check_constraints(spec, unit.members)
            if violations:
                unit.state = STATE_INCONSISTENT
                unit.reasons.extend(violations)
            else:
                unit.state = STATE_COMPLETE
        units.append(unit)
        if unit.state in (STATE_COMPLETE, STATE_PENDING) and unit.members:
            counting.append(list(unit.members))
        else:
            # Inconsistent / identity-ambiguous links: both ends stay
            # visible as independent fragment-level units.
            for member in unit.members:
                counting.append([member])

    linked_ids = set(assigned)
    for view in views:
        if view.fragment_id not in linked_ids:
            units.append(LinkUnit(
                unit_id=f"unit-{view.fragment_id}",
                link_id=None,
                molecule_id=None,
                state=STATE_SINGLETON,
                members=[view],
            ))
            counting.append([view])
    return units, counting


def combination_resolution(
    counting_units: list[list[FragmentView]],
    groups: list[set[str]],
) -> tuple[dict[int, dict[frozenset, str]], list[dict]]:
    """Combine member matches inside each multi-fragment unit.

    Within an ambiguity group, the interpretation that survives is the
    intersection of the candidates each member supports: a single shared
    candidate resolves the unit; disjoint member matches are a
    molecule-level conflict and support nothing.
    """
    resolved: dict[int, dict[frozenset, str]] = {}
    conflicts: list[dict] = []
    for index, members in enumerate(counting_units):
        if len(members) < 2:
            continue
        for group in groups:
            member_sets = [
                set(m.alt_variant_ids) & group for m in members
            ]
            non_empty = [s for s in member_sets if s]
            if len(non_empty) < 2:
                continue
            common = set.intersection(*non_empty)
            if len(common) == 1:
                resolved.setdefault(index, {})[frozenset(group)] = next(
                    iter(common))
            elif not common:
                conflicts.append({
                    "unit_index": index,
                    "group": sorted(group),
                    "candidate_variant_ids": sorted(
                        set.union(*non_empty)),
                    "reason": "members_support_disjoint_candidates",
                    "detail": "linked members support mutually exclusive "
                              "candidates; the molecule abstains",
                })
    return resolved, conflicts


def molecule_support_counts(
    counting_units: list[list[FragmentView]],
    unit_ids: list[str],
    resolved: dict[int, dict[frozenset, str]],
    group_of: dict[str, frozenset],
    conflicts: list[dict] | None = None,
) -> dict[str, dict]:
    """Per-candidate molecule-level support, traceable to unit members."""
    counts: dict[str, dict] = {}

    def entry(vid: str) -> dict:
        return counts.setdefault(vid, {
            "alt_molecules": 0,
            "independent_alt_molecules": 0,
            "shared_ambiguous_molecules": 0,
            "linkage_resolved_molecules": 0,
            "molecule_ids": [],
        })

    abstaining: dict[int, list[frozenset]] = {}
    for conflict in conflicts or []:
        abstaining.setdefault(conflict["unit_index"], []).append(
            frozenset(conflict["group"]))

    for index, members in enumerate(counting_units):
        support: set[str] = set()
        for member in members:
            support.update(member.alt_variant_ids)
        if not support:
            continue
        for group in abstaining.get(index, []):
            # Molecule-level conflict: the unit abstains from every
            # candidate in the affected group.
            support -= group
        if not support:
            continue
        unit_resolution = resolved.get(index, {})
        for vid in sorted(support):
            group = group_of.get(vid)
            rec = entry(vid)
            rec["alt_molecules"] += 1
            rec["molecule_ids"].append(unit_ids[index])
            if group is not None and group in unit_resolution:
                if unit_resolution[group] == vid:
                    rec["independent_alt_molecules"] += 1
                    rec["linkage_resolved_molecules"] += 1
                # Resolved to another candidate: no support here.
                continue
            if group is not None and len(support & group) > 1:
                rec["shared_ambiguous_molecules"] += 1
            else:
                rec["independent_alt_molecules"] += 1
    return counts


def counting_unit_ids(
    units: list[LinkUnit],
    counting_units: list[list[FragmentView]],
) -> list[str]:
    """Stable unit id for each counting unit (linked unit or singleton)."""
    ids: list[str] = []
    for members in counting_units:
        owner = next(
            (
                unit
                for unit in units
                if unit.state != STATE_SINGLETON
                and any(m.fragment_id == members[0].fragment_id
                        for m in unit.members)
            ),
            None,
        )
        ids.append(
            owner.unit_id if owner is not None
            else f"unit-{members[0].fragment_id}"
        )
    return ids


def recount_payload(payload: dict, specs: list[LinkSpec]) -> dict:
    """Re-resolve linkage under (possibly alternative) link declarations
    using the per-fragment detail stored in a published payload.

    Lets a researcher compare how another linkage rule would change
    molecule-level support without re-running the analysis and without
    touching the immutable base version.
    """
    views = [
        FragmentView.from_detail(detail)
        for detail in payload.get("fragment_linkage_detail", [])
    ]
    units, counting_units = resolve_links(specs, views)
    unit_ids = counting_unit_ids(units, counting_units)
    groups = [
        set(group["candidate_variant_ids"])
        for group in payload.get("ambiguity_groups", [])
    ]
    resolved, conflicts = combination_resolution(counting_units, groups)
    group_of = {
        vid: frozenset(group)
        for group in groups
        for vid in group
    }
    counts = molecule_support_counts(counting_units, unit_ids, resolved,
                                     group_of)
    return {
        "units": [unit.to_dict() for unit in units],
        "molecule_counts": counts,
        "molecule_conflicts": conflicts,
    }
