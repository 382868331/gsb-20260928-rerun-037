from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .models import Observation, Parameters

_CONSUMES_READ = {"M", "I", "=", "X"}
_CONSUMES_REFERENCE = {"M", "D", "N", "=", "X"}
_SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400


@dataclass(frozen=True)
class Alignment:
    template_id: str
    chromosome: str
    reference_start0: int
    mapping_quality: int
    is_reverse: bool
    cigar: tuple[tuple[int, str], ...]
    sequence: str
    qualities: tuple[int, ...]


def parse_sam(path: str | Path, parameters: Parameters) -> list[Alignment]:
    alignments: list[Alignment] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                raise ValueError("SAM records must contain at least 11 tab-separated fields")
            flag = int(fields[1])
            if flag & _SKIP_FLAGS:
                continue
            if parameters.skip_duplicates and flag & 0x400:
                continue
            mapping_quality = int(fields[4])
            if mapping_quality < parameters.min_mapping_quality:
                continue
            alignments.append(
                Alignment(
                    template_id=fields[0],
                    chromosome=fields[2],
                    reference_start0=int(fields[3]) - 1,
                    mapping_quality=mapping_quality,
                    is_reverse=bool(flag & 0x10),
                    cigar=tuple(_parse_cigar(fields[5])),
                    sequence=fields[9].upper(),
                    qualities=tuple(ord(char) - 33 for char in fields[10]),
                )
            )
    _validate_alignments(alignments)
    return alignments


def _validate_alignments(alignments: list[Alignment]) -> None:
    for alignment in alignments:
        ref0 = alignment.reference_start0
        read0 = 0
        for length, operation in alignment.cigar:
            if operation in _CONSUMES_READ:
                if read0 + length > len(alignment.sequence) or read0 + length > len(
                    alignment.qualities
                ):
                    raise ValueError(
                        f"CIGAR extends beyond SEQ/QUAL for read {alignment.template_id}"
                    )
            if operation in _CONSUMES_REFERENCE:
                ref0 += length


def _validate_reference_bounds(alignments: list[Alignment], reference_length: int) -> None:
    for alignment in alignments:
        if alignment.reference_start0 < 0:
            raise ValueError(f"Read {alignment.template_id} starts before the reference")
        ref0 = alignment.reference_start0
        for length, operation in alignment.cigar:
            if operation in _CONSUMES_REFERENCE:
                ref0 += length
        if ref0 > reference_length:
            raise ValueError(f"Read {alignment.template_id} extends beyond the reference")


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    if cigar == "*":
        raise ValueError("Aligned read without CIGAR cannot be evaluated")
    operations: list[tuple[int, str]] = []
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
        else:
            if not number or char not in "MIDNSHP=X":
                raise ValueError(f"Unsupported CIGAR operation in: {cigar}")
            operations.append((int(number), char))
            number = ""
    if number:
        raise ValueError(f"Truncated CIGAR: {cigar}")
    return operations


def observations_and_coverage(
    alignments: list[Alignment], reference: str
) -> tuple[list[Observation], dict[int, dict[str, int]], list[tuple[str, int, int]]]:
    _validate_reference_bounds(alignments, len(reference))
    observations: list[Observation] = []
    base_quality_by_position: dict[int, dict[str, int]] = {}
    orientation_by_position: dict[int, dict[str, tuple[int, bool]]] = {}
    intervals: list[tuple[str, int, int]] = []

    for alignment in alignments:
        ref0 = alignment.reference_start0
        read0 = 0
        covered_start = ref0
        chromosome = alignment.chromosome

        for length, operation in alignment.cigar:
            if operation in {"P", "H"}:
                continue
            if operation == "S":
                read0 += length
                continue
            if operation == "D":
                if ref0 < 1 or ref0 + length > len(reference):
                    raise ValueError("Deletion is outside reference bounds")
                anchor = reference[ref0 - 1]
                anchor_quality = (
                    alignment.qualities[read0 - 1] if read0 > 0 else 0
                )
                observations.append(
                    Observation(
                        template_id=alignment.template_id,
                        chromosome=chromosome,
                        position=ref0,
                        ref=anchor + reference[ref0 : ref0 + length],
                        alt=anchor,
                        variant_type="deletion",
                        is_reverse=alignment.is_reverse,
                        base_quality=anchor_quality,
                        mapping_quality=alignment.mapping_quality,
                    )
                )
                for deleted_ref0 in range(ref0, ref0 + length):
                    _record_base_quality(
                        base_quality_by_position,
                        orientation_by_position,
                        deleted_ref0,
                        alignment.template_id,
                        anchor_quality,
                        alignment.is_reverse,
                    )
                ref0 += length
                continue
            if operation == "N":
                ref0 += length
                continue

            if operation in _CONSUMES_REFERENCE:
                for _ in range(length):
                    _record_base_quality(
                        base_quality_by_position,
                        orientation_by_position,
                        ref0,
                        alignment.template_id,
                        alignment.qualities[read0],
                        alignment.is_reverse,
                    )
                    if operation in {"M", "=", "X"}:
                        read_base = alignment.sequence[read0]
                        ref_base = reference[ref0]
                        if read_base in "ACGT" and ref_base in "ACGT" and read_base != ref_base:
                            observations.append(
                                Observation(
                                    template_id=alignment.template_id,
                                    chromosome=chromosome,
                                    position=ref0 + 1,
                                    ref=ref_base,
                                    alt=read_base,
                                    variant_type="substitution",
                                    is_reverse=alignment.is_reverse,
                                    base_quality=alignment.qualities[read0],
                                    mapping_quality=alignment.mapping_quality,
                                )
                            )
                    if operation in _CONSUMES_READ:
                        read0 += 1
                    ref0 += 1
            elif operation == "I":
                if ref0 <= 0 or ref0 > len(reference):
                    raise ValueError("Insertion is outside reference bounds")
                anchor = reference[ref0 - 1]
                inserted = alignment.sequence[read0 : read0 + length]
                qualities = alignment.qualities[read0 : read0 + length]
                boundary = alignment.qualities[read0 - 1 : read0] + qualities
                observations.append(
                    Observation(
                        template_id=alignment.template_id,
                        chromosome=chromosome,
                        position=ref0,
                        ref=anchor,
                        alt=anchor + inserted,
                        variant_type="insertion",
                        is_reverse=alignment.is_reverse,
                        base_quality=min(boundary),
                        mapping_quality=alignment.mapping_quality,
                    )
                )
                read0 += length
        intervals.append((alignment.template_id, covered_start, ref0))

    return observations, base_quality_by_position, orientation_by_position, intervals


def _record_base_quality(
    quality_map: dict[int, dict[str, int]],
    orientation_map: dict[int, dict[str, tuple[int, bool]]],
    ref0: int,
    template_id: str,
    quality: int,
    is_reverse: bool,
) -> None:
    current = quality_map.setdefault(ref0, {})
    current[template_id] = max(current.get(template_id, -1), quality)
    orientations = orientation_map.setdefault(ref0, {})
    if quality > orientations.get(template_id, (-1, False))[0]:
        orientations[template_id] = (quality, is_reverse)
