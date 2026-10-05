"""Aggregate per-fragment evidence into candidate variants.

One DNA fragment contributes at most one vote at an anchor position, even
when both mates overlap it. Low-quality alt observations stay visible (they
lower confidence) but never count as passing support.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import Parameters
from .evidence import CandidateEvidence, fisher_two_sided, score_candidate
from .fasta import Reference
from .normalize import Variant
from .samio import Fragment


@dataclass
class CoveragePosition:
    chrom: str
    pos: int
    fragment_depth: int
    passing_depth: int


@dataclass
class CallResult:
    candidates: list[CandidateEvidence]
    coverage: list[CoveragePosition]
    low_coverage: list[tuple[str, int, int, int]]  # chrom,start,end,min_depth
    ambiguity_groups: list[dict] = field(default_factory=list)
    fragment_sources: list[dict] = field(default_factory=list)
    evidence_units: list[dict] = field(default_factory=list)
    combination_search: dict | None = None


def _overlaps(a: Variant, b: Variant) -> bool:
    return a.chrom == b.chrom and a.pos < b.pos + len(b.ref) and b.pos < a.pos + len(a.ref)


def _spans(frag: Fragment, chrom: str, anchor: int) -> bool:
    if chrom not in frag.chroms and len(frag.chroms) > 1:
        return False
    return any(start <= anchor <= end for start, end in frag.covered)


def _fragment_state(frag: Fragment, key: tuple, anchor: int) -> str:
    if key in frag.alt:
        return "alt_low" if frag.alt[key].low_quality else "alt_pass"
    if anchor in frag.conflict_positions:
        return "conflict"
    for other_key, obs in frag.alt.items():
        if other_key == key:
            continue
        variant = obs.variant
        if variant.pos <= anchor < variant.pos + len(variant.ref):
            return "other_alt"
    if anchor in frag.ref_positions:
        return "ref"
    return "uninformative"


def _depth_arrays(
    fragments: list[Fragment], references: dict[str, Reference],
    fragment_units: dict[str, str] | None = None,
) -> dict[str, tuple[list[int], list[int]]]:
    fragment_units = fragment_units or {}
    raw = {
        name: [0] * (len(ref.sequence) + 2)
        for name, ref in references.items()
    }
    counted_intervals: set[tuple[str, str, int, int]] = set()
    for frag in fragments:
        chroms = frag.chroms
        unit_key = fragment_units.get(frag.evidence_id, frag.evidence_id)
        for start, end in frag.covered:
            chrom = next(iter(chroms)) if len(chroms) == 1 else None
            if chrom is None:
                # Multi-contig fragment: spans cannot be assigned safely.
                continue
            interval_key = (unit_key, chrom, start, end)
            if interval_key in counted_intervals:
                continue
            counted_intervals.add(interval_key)
            raw[chrom][start] += 1
            raw[chrom][end + 1] -= 1
    out: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, diff in raw.items():
        depth = [0] * len(diff)
        running = 0
        for i in range(1, len(diff) - 1):
            running += diff[i]
            depth[i] = running
        out[chrom] = (depth, depth)  # MAPQ-filtered fragments already excluded
    return out


def call_variants(
    fragments: list[Fragment],
    references: dict[str, Reference],
    params: Parameters,
    association_relations: list | None = None,
    combination_constraints: list | None = None,
) -> CallResult:
    fragment_sources = [
        {
            "fragment_id": frag.evidence_id,
            "qname": frag.qname,
            "source_id": frag.source_id,
            "source_name": frag.source_name or frag.source_id,
            "batch_id": frag.batch_id or frag.source_id,
        }
        for frag in fragments
    ]
    grouped: dict[tuple, list] = {}
    variant_lookup: dict[tuple, Variant] = {}
    for frag in fragments:
        for key, obs in frag.alt.items():
            grouped.setdefault(key, []).append(obs)
            variant_lookup[key] = obs.variant

    candidate_results: list[CandidateEvidence] = []
    resolved_units = []
    if association_relations:
        from .association import (
            MOLECULE_CONFLICT,
            RESOLVED,
            reconcile_associations,
        )
        resolved_units = [
            unit.to_dict()
            for unit in reconcile_associations(fragments, association_relations)
            if unit.status in (RESOLVED, MOLECULE_CONFLICT)
        ]

    def association_context(fragment_id: str, variant_id: str):
        for unit in resolved_units:
            if fragment_id not in unit["fragment_ids"]:
                continue
            if variant_id in unit.get("conflicting_candidates", []):
                return "conflict", unit["association_id"]
            if variant_id in unit.get("excluded_candidates", []):
                return "excluded", None
            if variant_id in unit.get("linked_candidates", []):
                return "linked", unit["association_id"]
        return "unlinked", None

    fragment_units = {
        fragment_id: unit["association_id"]
        for unit in resolved_units
        for fragment_id in unit["fragment_ids"]
    }
    excluded_by_unit = {
        unit["association_id"]: set(unit.get("excluded_candidates", []))
        for unit in resolved_units
    }
    fragment_excluded: dict[str, set[str]] = {}
    for unit in resolved_units:
        excluded = set(unit.get("excluded_candidates", []))
        for fragment_id, candidates in unit.get("single_fragment_candidates", {}).items():
            fragment_excluded.setdefault(fragment_id, set()).update(
                excluded & set(candidates))

    for key, observations in grouped.items():
        variant = variant_lookup[key]
        anchor = variant.pos
        chrom = variant.chrom

        states_by_unit: dict[str, dict] = {}
        n_alt_pass = n_alt_low = n_ref = n_other = n_conflict = 0
        n_independent = n_shared = 0
        fwd_alt = rev_alt = fwd_ref = rev_ref = 0
        depth = passing_depth = 0
        counted_depth_keys: set[str] = set()
        quals: list[int] = []
        fragment_ids: list[str] = []
        shared_with: dict[str, list[str]] = {}

        for frag in fragments:
            if not _spans(frag, chrom, anchor):
                continue
            depth_key = fragment_units.get(frag.evidence_id, frag.evidence_id)
            if depth_key not in counted_depth_keys:
                counted_depth_keys.add(depth_key)
                depth += 1
            association_state, unit_id = association_context(
                frag.evidence_id, variant.vid)
            support_key = unit_id or frag.evidence_id
            if association_state == "conflict":
                state = "conflict"
            elif association_state == "excluded":
                state = "other_alt"
            else:
                state = _fragment_state(frag, key, anchor)
            if state == "other_alt":
                excluded_here = {
                    observation.variant.vid
                    for other_key, observation in frag.alt.items()
                    if other_key != key and _overlaps(variant, observation.variant)
                } & fragment_excluded.get(frag.evidence_id, set())
                if excluded_here:
                    state = "uninformative"
            if state == "uninformative":
                continue
            entry = states_by_unit.setdefault(support_key, {
                "state": state, "fragment": frag, "unit_id": unit_id,
                "has_alt": False})
            # Merge endpoint states within one molecule:
            # - conflict (same-site contradictory alleles) is decisive and
            #   forces the molecule to abstain at this site;
            # - any endpoint with ALT at this compatible site means the
            #   molecule supports it once, regardless of other-site conflicts.
            if entry["state"] != "conflict":
                if state in ("alt_pass", "alt_low"):
                    entry.update(state=state, fragment=frag, unit_id=unit_id,
                                 has_alt=True)
                elif not entry["has_alt"]:
                    rank = {"other_alt": 3, "alt_low": 2, "alt_pass": 1, "ref": 0}
                    if rank[state] > rank[entry["state"]]:
                        entry.update(state=state, fragment=frag, unit_id=unit_id)

        passing_depth = len(states_by_unit)
        for support_key, entry in states_by_unit.items():
            state, frag, unit_id = (
                entry["state"], entry["fragment"], entry["unit_id"])
            if state in ("alt_pass", "alt_low"):
                fragment_ids.append(support_key)
            if state == "alt_pass":
                n_alt_pass += 1
                other_ids = [
                    other.variant.vid
                    for other_key, other in frag.alt.items()
                    if other_key != key and _overlaps(variant, other.variant)
                    and other.variant.vid not in excluded_by_unit.get(unit_id, set())
                ]
                if other_ids:
                    n_shared += 1
                    shared_with[support_key] = sorted(set(other_ids))
                else:
                    n_independent += 1
                quals.extend(frag.alt[key].base_quals)
                fwd_alt += frag.strand != "R"
                rev_alt += frag.strand == "R"
            elif state == "alt_low":
                n_alt_low += 1
                other_ids = [
                    other.variant.vid
                    for other_key, other in frag.alt.items()
                    if other_key != key and _overlaps(variant, other.variant)
                    and other.variant.vid not in excluded_by_unit.get(unit_id, set())
                ]
                if other_ids:
                    n_shared += 1
                else:
                    n_independent += 1
            elif state == "ref":
                n_ref += 1
                fwd_ref += frag.strand != "R"
                rev_ref += frag.strand == "R"
            elif state == "other_alt":
                n_other += 1
            elif state == "conflict":
                n_conflict += 1

        informative = n_alt_pass + n_ref
        alt_fraction = n_alt_pass / informative if informative else 0.0
        mean_q = sum(quals) / len(quals) if quals else 0.0
        strand_p = fisher_two_sided(fwd_alt, rev_alt, fwd_ref, rev_ref)

        ev = CandidateEvidence(
            variant_id=variant.vid,
            chrom=chrom,
            pos=anchor,
            ref=variant.ref,
            alt=variant.alt,
            kind=variant.kind,
            alt_fragments=n_alt_pass + n_alt_low,
            alt_forward=fwd_alt,
            alt_reverse=rev_alt,
            alt_mean_base_quality=round(mean_q, 2),
            low_quality_alt_fragments=n_alt_low,
            alt_fragment_ids=fragment_ids,
            independent_alt_fragments=n_independent,
            shared_ambiguous_alt_fragments=n_shared,
            ref_fragments=n_ref,
            other_alt_fragments=n_other,
            conflicting_overlap_fragments=n_conflict,
            strand_bias_p=round(strand_p, 8),
            depth=depth,
            passing_depth=passing_depth,
            alt_fraction=round(alt_fraction, 6),
            informative_depth=informative,
        )
        score_candidate(
            ev,
            min_alt_count=params.min_alt_count,
            min_alt_fraction=params.min_alt_fraction,
            min_coverage=params.min_coverage_for_calling,
            strand_bias_p_threshold=params.strand_bias_p_threshold,
            low_quality_fraction_threshold=params.low_quality_fraction_threshold,
        )
        ev.shared_with = shared_with
        candidate_results.append(ev)

    candidate_results.sort(key=lambda e: (e.chrom, e.pos, e.alt))
    depth_maps = _depth_arrays(fragments, references, fragment_units)
    coverage, low = _coverage_report(depth_maps, params.low_coverage_depth)
    ambiguity_groups = _ambiguity_groups(candidate_results)
    evidence_units = []
    if association_relations:
        from .association import reconcile_associations
        evidence_units = [
            unit.to_dict()
            for unit in reconcile_associations(fragments, association_relations)
        ]
    combination_search = None
    if combination_constraints:
        from .combination import search_combinations
        combination_search = search_combinations(
            fragments, association_relations or [], combination_constraints
        ).to_dict()
    return CallResult(
        candidates=candidate_results,
        coverage=coverage,
        low_coverage=low,
        ambiguity_groups=ambiguity_groups,
        fragment_sources=fragment_sources,
        evidence_units=evidence_units,
        combination_search=combination_search,
    )


def _ambiguity_groups(candidates: list[CandidateEvidence]) -> list[dict]:
    adjacency = {ev.variant_id: set() for ev in candidates}
    fragment_members: dict[str, set[str]] = {}
    for ev in candidates:
        for fragment_id, others in getattr(ev, "shared_with", {}).items():
            adjacency.setdefault(ev.variant_id, set()).update(others)
            for other in others:
                adjacency.setdefault(other, set()).add(ev.variant_id)
            fragment_members.setdefault(fragment_id, {ev.variant_id}).add(ev.variant_id)
            fragment_members[fragment_id].update(others)

    seen: set[str] = set()
    groups: list[dict] = []
    for variant_id in sorted(adjacency):
        if variant_id in seen:
            continue
        stack = [variant_id]
        members: set[str] = set()
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            members.add(current)
            stack.extend(adjacency.get(current, set()) - seen)
        if len(members) < 2:
            continue
        records = []
        for fragment_id, linked in sorted(fragment_members.items()):
            if linked & members:
                records.append({
                    "fragment_id": fragment_id,
                    "supports": sorted(linked & members),
                })
        groups.append({
            "group_id": f"ambig-{len(groups) + 1:04d}",
            "candidate_variant_ids": sorted(members),
            "shared_fragments": records,
            "shared_fragment_count": len(records),
            "resolution": "ambiguous",
        })
    return groups


def _coverage_report(
    depth_maps: dict[str, tuple[list[int], list[int]]], cutoff: int
):
    coverage: list[CoveragePosition] = []
    low: list[tuple[str, int, int, int]] = []
    for chrom, (raw, passing) in sorted(depth_maps.items()):
        run_start: int | None = None
        min_depth = 0
        length = len(raw) - 2
        for pos in range(1, length + 1):
            d = raw[pos]
            coverage.append(
                CoveragePosition(
                    chrom=chrom,
                    pos=pos,
                    fragment_depth=d,
                    passing_depth=passing[pos],
                )
            )
            if d < cutoff:
                if run_start is None:
                    run_start = pos
                    min_depth = d
                else:
                    min_depth = min(min_depth, d)
            elif run_start is not None:
                low.append((chrom, run_start, pos - 1, min_depth))
                run_start = None
        if run_start is not None:
            low.append((chrom, run_start, length, min_depth))
    return coverage, low
