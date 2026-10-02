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
    fragments: list[Fragment], references: dict[str, Reference]
) -> dict[str, tuple[list[int], list[int]]]:
    raw = {
        name: [0] * (len(ref.sequence) + 2)
        for name, ref in references.items()
    }
    for frag in fragments:
        chroms = frag.chroms
        for start, end in frag.covered:
            chrom = next(iter(chroms)) if len(chroms) == 1 else None
            if chrom is None:
                # Multi-contig fragment: spans cannot be assigned safely.
                continue
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

    candidates: list[CandidateEvidence] = []
    for key, observations in grouped.items():
        variant = variant_lookup[key]
        anchor = variant.pos
        chrom = variant.chrom

        n_alt_pass = n_alt_low = n_ref = n_other = n_conflict = 0
        n_independent = n_shared = 0
        fwd_alt = rev_alt = fwd_ref = rev_ref = 0
        depth = passing_depth = 0
        quals: list[int] = []
        fragment_ids: list[str] = []
        shared_with: dict[str, list[str]] = {}

        for frag in fragments:
            if not _spans(frag, chrom, anchor):
                continue
            depth += 1
            state = _fragment_state(frag, key, anchor)
            if state == "uninformative":
                continue
            passing_depth += 1
            if state == "alt_pass":
                n_alt_pass += 1
                fragment_ids.append(frag.evidence_id)
                other_ids = [
                    other.variant.vid
                    for other_key, other in frag.alt.items()
                    if other_key != key and _overlaps(variant, other.variant)
                ]
                if other_ids:
                    n_shared += 1
                    shared_with[frag.evidence_id] = sorted(set(other_ids))
                else:
                    n_independent += 1
                quals.extend(frag.alt[key].base_quals)
                if frag.strand == "R":
                    rev_alt += 1
                else:
                    fwd_alt += 1
            elif state == "alt_low":
                n_alt_low += 1
                fragment_ids.append(frag.evidence_id)
                other_ids = [
                    other.variant.vid
                    for other_key, other in frag.alt.items()
                    if other_key != key and _overlaps(variant, other.variant)
                ]
                if other_ids:
                    n_shared += 1
                    shared_with[frag.evidence_id] = sorted(set(other_ids))
                else:
                    n_independent += 1
            elif state == "ref":
                n_ref += 1
                if frag.strand == "R":
                    rev_ref += 1
                else:
                    fwd_ref += 1
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
        candidates.append(ev)

    candidates.sort(key=lambda e: (e.chrom, e.pos, e.alt))
    depth_maps = _depth_arrays(fragments, references)
    coverage, low = _coverage_report(depth_maps, params.low_coverage_depth)
    ambiguity_groups = _ambiguity_groups(candidates)
    return CallResult(
        candidates=candidates,
        coverage=coverage,
        low_coverage=low,
        ambiguity_groups=ambiguity_groups,
        fragment_sources=fragment_sources,
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
