"""Allele normalization and stable variant identity.

All callers (this stage and future annotation) derive the same canonical key
from the parsimonious, left-aligned allele representation, so indels that are
spelled differently by different reads collapse to one identity.
"""

from __future__ import annotations

from dataclasses import dataclass

from .fasta import Reference


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int          # 1-based position of the first (anchor) base
    ref: str          # reference allele, uppercase
    alt: str          # alternate allele, uppercase
    kind: str         # "SNV" | "INS" | "DEL" | "MNP"

    @property
    def vid(self) -> str:
        if self.kind == "SNV":
            return f"{self.chrom}:g.{self.pos}{self.ref}>{self.alt}"
        if self.kind == "INS":
            return f"{self.chrom}:g.{self.pos}_{self.pos + 1}ins{self.alt[1:]}"
        if self.kind == "DEL":
            end = self.pos + len(self.ref) - 1
            return f"{self.chrom}:g.{self.pos}_{end}del{self.ref[1:]}"
        end = self.pos + len(self.ref) - 1
        return f"{self.chrom}:g.{self.pos}_{end}{self.ref}>{self.alt}"

    @property
    def canonical_key(self) -> tuple[str, int, str, str]:
        return (self.chrom, self.pos, self.ref, self.alt)


def classify_alleles(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if ref[0] == alt[0] and len(alt) > len(ref):
        return "INS"
    if ref[0] == alt[0] and len(ref) > len(alt):
        return "DEL"
    return "MNP"


def normalize(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    reference: Reference,
) -> Variant:
    """Trim redundant suffix/prefix bases then left-align repeats.

    ``pos`` is 1-based and points at ``ref[0]``. Insertions/deletions are
    expected to carry a single anchor base on each allele.
    """
    ref = ref.upper()
    alt = alt.upper()
    if not ref or not alt:
        raise ValueError("alleles must be non-empty and anchor-padded")

    def base_at(position_1based: int) -> str | None:
        idx = position_1based - 1
        if 0 <= idx < len(reference.sequence):
            return reference.sequence[idx]
        return None

    # 1. trim shared trailing bases while each allele keeps its anchor + at
    #    least one event base (the single-anchor case is handled by left-align).
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    # 2. left-align pure indels across identical preceding bases (tandems).
    if len(ref) == 1 or len(alt) == 1:
        while True:
            moved = ref[-1] if len(ref) > 1 else alt[-1]
            preceding = base_at(pos - 1)
            if preceding is None or preceding != moved:
                break
            if len(ref) > 1 and len(alt) > 1:
                break
            ref = preceding + ref[:-1]
            alt = preceding + alt[:-1]
            pos -= 1

    # 3. trim shared leading bases (MNP-style parsimony) while keeping anchor.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    kind = classify_alleles(ref, alt)
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, kind=kind)
