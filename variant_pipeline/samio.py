"""Minimal SAM alignment reader (SAM text, no external dependencies).

Orientation convention:

Every accepted record is interpreted on the forward reference axis. POS is
the 1-based reference position of the first consumed reference base and the
CIGAR (M/=/X, I, D/N, S/H/P) is walked in forward reference order. FLAG 0x10
only labels the strand for strand-bias statistics and never changes base
coordinates or quality order. SEQ/QUAL are therefore compared to the
reference directly and are never reverse complemented again: a fully
matching reverse read stored in reference orientation contributes zero
variants, the same as a forward read.

Clips (S/H/P) produce no evidence. Secondary, supplementary and unmapped
records are ignored. Mate overlap is merged per DNA fragment so reads are
never double counted; conflicting overlap alleles abstain at that position.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from .fasta import Reference
from .normalize import Variant, normalize

FLAG_REVERSE = 0x10
FLAG_FIRST = 0x40
FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


@dataclass
class Observation:
    """One fragment's allele call at one canonical event."""

    variant: Variant
    base_quals: tuple[int, ...]
    strand: str            # fragment strand: 'F' forward / 'R' reverse
    low_quality: bool      # any supporting base below the quality gate


@dataclass
class Fragment:
    qname: str
    strand: str
    source_id: str = "default"
    source_name: str | None = None
    batch_id: str | None = None
    chroms: set[str] = field(default_factory=set)
    alt: dict[tuple, Observation] = field(default_factory=dict)
    ref_positions: set[int] = field(default_factory=set)
    conflict_positions: set[int] = field(default_factory=set)
    covered: list[tuple[int, int]] = field(default_factory=list)
    mapq_max: int = 0

    @property
    def evidence_id(self) -> str:
        """Stable identity that keeps same qname from different sources apart."""
        return f"{self.source_id}::{self.qname}"


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    if cigar in ("", "*"):
        return []
    parts = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if "".join(n + op for n, op in parts) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar!r}")
    return [(int(n), op) for n, op in parts]


@dataclass(frozen=True)
class AlignmentRecord:
    qname: str
    flag: int
    rname: str
    pos: int            # 1-based
    mapq: int
    seq: str
    quals: tuple[int, ...]
    cigar: str

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & FLAG_REVERSE)


def _parse_record(columns: list[str]) -> AlignmentRecord:
    if len(columns) < 11:
        raise ValueError("SAM record has fewer than 11 mandatory fields")
    qname, flag, rname, pos, mapq, cigar, _rnext, _pnext, _tlen, seq, qual = columns[:11]
    if seq == "*":
        bases = ""
        quals = ()
    elif qual == "*":
        bases = seq.upper()
        quals = tuple([0] * len(bases))
    else:
        bases = seq.upper()
        quals = tuple(ord(c) - 33 for c in qual)
        if len(quals) != len(bases):
            raise ValueError(f"SEQ/QUAL length mismatch for {qname}")
    return AlignmentRecord(
        qname=qname,
        flag=int(flag),
        rname=rname,
        pos=int(pos),
        mapq=int(mapq),
        seq=bases,
        quals=quals,
        cigar=cigar,
    )


def read_sam(path: str | Path) -> tuple[dict[str, int], list[AlignmentRecord]]:
    """Return (@SQ length hints, usable primary alignment records)."""
    lengths: dict[str, int] = {}
    records: list[AlignmentRecord] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    name = length = None
                    for field_ in line.split("\t")[1:]:
                        if field_.startswith("SN:"):
                            name = field_[3:]
                        elif field_.startswith("LN:"):
                            length = int(field_[3:])
                    if name and length:
                        lengths[name] = length
                continue
            rec = _parse_record(line.split("\t"))
            if rec.flag & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY):
                continue
            if rec.rname == "*" or rec.pos <= 0 or not rec.seq:
                continue
            records.append(rec)
    return lengths, records


@dataclass
class _MateCall:
    strand: str
    span: tuple[int, int] | None
    alt_ranges: list[tuple[int, int, tuple, tuple[int, ...]]]
    variants: dict[tuple, Variant]
    ref_positions: set[int]


def _mate_call(rec: AlignmentRecord, reference: Reference) -> _MateCall:
    # SEQ/QUAL and CIGAR are already expressed on the forward reference axis
    # for both read orientations; only strand labelling uses FLAG 0x10.
    seq = rec.seq
    quals = list(rec.quals)
    qi = 0
    ri = rec.pos
    span_start = rec.pos
    span_end = rec.pos - 1
    alt_ranges: list[tuple[int, int, tuple, tuple[int, ...]]] = []
    variants: dict[tuple, Variant] = {}
    ref_positions: set[int] = set()
    last_match: tuple[int, int] | None = None  # (ref_pos, query_index)

    def add_alt(variant: Variant, bquals: tuple[int, ...]) -> None:
        end = variant.pos + len(variant.ref) - 1
        key = variant.canonical_key
        alt_ranges.append((variant.pos, end, key, bquals))
        variants[key] = variant

    for length, op in parse_cigar(rec.cigar):
        if op in ("M", "=", "X"):
            for _ in range(length):
                qb = seq[qi]
                rb = (
                    reference.sequence[ri - 1]
                    if 0 <= ri - 1 < len(reference.sequence)
                    else "N"
                )
                if qb in "ACGT" and rb in "ACGT":
                    if qb == rb:
                        ref_positions.add(ri)
                    else:
                        add_alt(normalize(rec.rname, ri, rb, qb, reference), (quals[qi],))
                last_match = (ri, qi)
                qi += 1
                ri += 1
                span_end = ri - 1
        elif op == "I":
            if last_match is not None:
                anchor_ref_pos, anchor_qi = last_match
                anchor = seq[anchor_qi]
                inserted = seq[qi:qi + length]
                bquals = tuple([quals[anchor_qi]] + quals[qi:qi + length])
                add_alt(
                    normalize(
                        rec.rname,
                        anchor_ref_pos,
                        anchor,
                        anchor + inserted,
                        reference,
                    ),
                    bquals,
                )
            qi += length
        elif op in ("D", "N"):
            deleted = reference.sequence[ri - 1:ri - 1 + length]
            if last_match is not None:
                anchor_ref_pos, anchor_qi = last_match
                anchor = seq[anchor_qi]
                variant = normalize(
                    rec.rname,
                    anchor_ref_pos,
                    anchor + deleted,
                    anchor,
                    reference,
                )
                add_alt(variant, (quals[anchor_qi],))
            else:  # deletion starts at alignment start; anchor on following base
                anchor_pos = ri + length
                anchor_qi = qi
                anchor = seq[anchor_qi]
                variant = normalize(
                    rec.rname,
                    anchor_pos,
                    deleted + anchor,
                    anchor,
                    reference,
                )
                add_alt(variant, (quals[anchor_qi],))
            ri += length
            span_end = ri - 1
        elif op == "S":
            qi += length
        elif op in ("H", "P"):
            pass

    span = (span_start, span_end) if span_end >= span_start else None
    return _MateCall(
        strand="R" if rec.is_reverse else "F",
        span=span,
        alt_ranges=alt_ranges,
        variants=variants,
        ref_positions=ref_positions,
    )


def _mean(values) -> float:
    return sum(values) / len(values) if values else 0.0


def _allele_state_at(call: _MateCall, position: int):
    for start, end, key, _bquals in call.alt_ranges:
        if start <= position <= end:
            return key
    if position in call.ref_positions:
        return "REF"
    return None


def _allele_conflicts(a: _MateCall, b: _MateCall) -> set[int]:
    if a.span is None or b.span is None:
        return set()
    lo = max(a.span[0], b.span[0])
    hi = min(a.span[1], b.span[1])
    if lo > hi:
        return set()
    conflicts: set[int] = set()
    for position in range(lo, hi + 1):
        sa = _allele_state_at(a, position)
        sb = _allele_state_at(b, position)
        if sa is not None and sb is not None and sa != sb:
            conflicts.add(position)
    return conflicts


def build_fragments(
    records: list[AlignmentRecord],
    references: dict[str, Reference],
    min_base_quality: int,
    min_mapping_quality: int,
    source_id: str = "default",
    source_name: str | None = None,
    batch_id: str | None = None,
) -> list[Fragment]:
    grouped: dict[tuple[str, str], list[AlignmentRecord]] = {}
    for rec in records:
        grouped.setdefault((source_id, rec.qname), []).append(rec)

    fragments: list[Fragment] = []
    for (_record_source, qname), mates in grouped.items():
        usable = [m for m in mates if m.mapq >= min_mapping_quality]
        if not usable:
            continue
        read1 = next((m for m in usable if m.flag & FLAG_FIRST), usable[0])
        frag = Fragment(
            qname=qname,
            strand="R" if read1.is_reverse else "F",
            source_id=source_id,
            source_name=source_name,
            batch_id=batch_id,
            chroms={m.rname for m in usable if m.rname in references},)
        frag.mapq_max = max(m.mapq for m in usable)

        calls: list[_MateCall] = []
        for mate in usable:
            reference = references.get(mate.rname)
            if reference is None:
                continue
            call = _mate_call(mate, reference)
            if call.span is not None:
                frag.covered.append(call.span)
            calls.append(call)

        conflict: set[int] = set()
        for i in range(len(calls)):
            for j in range(i + 1, len(calls)):
                conflict |= _allele_conflicts(calls[i], calls[j])
        frag.conflict_positions = conflict

        merged: dict[tuple, Observation] = {}
        for call in calls:
            for start, end, key, bquals in call.alt_ranges:
                if any(p in conflict for p in range(start, end + 1)):
                    continue
                low_q = any(q < min_base_quality for q in bquals)
                existing = merged.get(key)
                if existing is None:
                    merged[key] = Observation(
                        variant=call.variants[key],
                        base_quals=bquals,
                        strand=call.strand,
                        low_quality=low_q,
                    )
                else:
                    if _mean(bquals) > _mean(existing.base_quals):
                        existing.base_quals = bquals
                    existing.low_quality = existing.low_quality and low_q
        frag.alt = merged

        ref_positions: set[int] = set()
        for call in calls:
            ref_positions |= call.ref_positions
        for obs in merged.values():
            for p in range(obs.variant.pos, obs.variant.pos + len(obs.variant.ref)):
                ref_positions.discard(p)
        ref_positions -= conflict
        frag.ref_positions = ref_positions
        fragments.append(frag)

    return fragments
